# -*- coding: utf-8 -*-
"""
VARIABLE STAR SEARCHER - MOTORE DI CALCOLO
==========================================
© 2026 Edoardo Battistel

Cerca stelle variabili NON ancora catalogate nelle curve di luce di TESS, Kepler e K2:
  1. stessi elenchi di stelle e stessi download di Exoplanet Searcher (motore_exo);
  2. periodogramma di Lomb-Scargle, controllo del periodo doppio (binarie a eclisse),
     forma della curva ripiegata -> tipo VSX proposto (EA, EB, EW, RRAB, RRC, DSCT,
     HADS, GDOR, ROT, BY, UV…);
  3. brillamenti (flare) nelle curve a cadenza breve;
  4. confronto con il catalogo AAVSO VSX (copia HEASARC): si tengono solo le sconosciute;
  5. per ogni variabile nuova: grafici, scheda PDF in inglese e i campi del modulo
     "New star" di VSX (vsx_submission.txt) + una riga in vsx_da_inviare.csv;
  6. pagina risultati.html con tutte le variabili nuove (si apre dal PC).
"""

import csv
import html
import json
import math
import re
import shutil
import time
import urllib.parse
import warnings
from datetime import datetime
from pathlib import Path

import threading

import numpy as np

import motore_exo as M

NOME_PROGRAMMA = "Variable Star Searcher"
VERSIONE = "1.2"
AUTORE = "Edoardo Battistel"
COPYRIGHT = f"© 2026 {AUTORE}"
M.INTESTAZIONI_HTTP = {"User-Agent": "VariableStarSearcher/1.2 (amateur astronomy software)"}

URL_HEASARC = "https://heasarc.gsfc.nasa.gov/xamin/vo/tap/sync"
RAGGIO_VSX_ARCSEC = 40.0          # un pixel TESS è 21": entro 2 pixel è "la stessa stella"

# --- soglie
FAP_MAX = 1e-8                    # probabilità di falso allarme del periodogramma
AMPIEZZA_MIN_MAG = 0.005          # 5 millimagnitudini: sotto è troppo facile confondersi con i difetti
SNR_MIN = 12.0                    # ampiezza / rumore della curva ripiegata
CROWD_MIN = 0.7                   # meno del 70 % della luce dalla stella: troppo contaminata
PERIODI_SISTEMATICI = (13.7, 6.85, 27.4, 1.0)   # orbita di TESS e giorno: periodi sospetti

_scrittura = threading.Lock()      # più stelle analizzate insieme scrivono sugli stessi file
FILE_RIASSUNTO = "riassunto_stelle.csv"
FILE_CANDIDATI = "variabili_nuove.csv"
FILE_VSX = "vsx_da_inviare.csv"
FILE_PAGINA = "risultati.html"
CARTELLA_CANDIDATI = "variabili"

OFFSET_BJD = {"tess2": 2457000.0, "tessffi": 2457000.0, "kepler": 2454833.0, "k2": 2454833.0}
BANDA_VSX = {"tess2": "TESS", "tessffi": "TESS", "kepler": "Kp", "k2": "Kp"}
NOMI_TIPI = {  # tipo VSX -> descrizione in inglese per la scheda
    "EA": "Algol-type eclipsing binary", "EB": "Beta Lyrae-type eclipsing binary",
    "EW": "W UMa-type eclipsing binary", "ELL": "ellipsoidal variable", "RRAB": "RR Lyrae (fundamental mode)",
    "RRC": "RR Lyrae (first overtone)", "DSCT": "Delta Scuti", "HADS": "high-amplitude Delta Scuti",
    "GDOR": "Gamma Doradus", "ROT": "rotational (spotted) variable", "BY": "BY Draconis-type (spotted dwarf)",
    "UV": "UV Ceti-type flare star", "VAR": "variable of unspecified type",
}


# =============================================================================
# CURVA DI LUCE
# =============================================================================
def leggi(percorsi):
    """Curva (tempo, flusso normalizzato, errore) e informazioni, con contaminazione e cadenza."""
    t, f, e, info = M.leggi_curva(percorsi)
    try:
        from astropy.io import fits
        with fits.open(percorsi[0], memmap=False) as h:
            for k in ("CROWDSAP",):
                v = h[1].header.get(k)
                if isinstance(v, (int, float)) and math.isfinite(float(v)):
                    info["crowdsap"] = float(v)
            for k in ("SECTOR", "QUARTER", "CAMPAIGN", "OBJECT", "TICID", "KEPLERID"):
                v = h[0].header.get(k)
                if v is not None:
                    info[k.lower()] = v
    except Exception:
        pass
    dt = np.diff(t)
    info["cadenza"] = float(np.median(dt[dt > 0])) if len(dt) else 0.02
    return t, f, e, info


def rumore_punti(f):
    """Rumore da punto a punto (non risente della variabilità lenta)."""
    d = np.diff(f)
    return 1.4826 * np.median(np.abs(d - np.median(d))) / math.sqrt(2) if len(d) > 10 else float(np.std(f))


# =============================================================================
# PERIODOGRAMMA
# =============================================================================
def _ls_numpy(t, y, freq):
    """Lomb-Scargle semplice (riserva se manca astropy). Potenza normalizzata 0-1."""
    y = y - y.mean()
    var = np.sum(y * y)
    p = np.zeros(len(freq))
    for i, fr in enumerate(freq):
        w = 2 * np.pi * fr
        tau = math.atan2(np.sum(np.sin(2 * w * t)), np.sum(np.cos(2 * w * t))) / (2 * w)
        c, s = np.cos(w * (t - tau)), np.sin(w * (t - tau))
        p[i] = (np.sum(y * c) ** 2 / max(np.sum(c * c), 1e-12) + np.sum(y * s) ** 2 / max(np.sum(s * s), 1e-12)) / (2 * var)
    return p


def periodogramma(t, y, pmin, pmax):
    """(frequenze, potenze, periodo migliore, FAP)."""
    base = t.max() - t.min()
    try:
        from astropy.timeseries import LombScargle
        ls = LombScargle(t, y)
        freq, pot = ls.autopower(minimum_frequency=1 / pmax, maximum_frequency=1 / pmin, samples_per_peak=10)
        j = int(np.argmax(pot))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fap = float(ls.false_alarm_probability(pot[j], minimum_frequency=1 / pmax, maximum_frequency=1 / pmin))
    except ImportError:
        n = int(min(20000, 10 * base * (1 / pmin - 1 / pmax)))
        freq = np.linspace(1 / pmax, 1 / pmin, max(n, 500))
        pot = _ls_numpy(t, y, freq)
        j = int(np.argmax(pot))
        nf = base * (1 / pmin - 1 / pmax)
        fap = float(min(1.0, nf * (1 - pot[j]) ** ((len(t) - 3) / 2)))
    return np.asarray(freq), np.asarray(pot), float(1 / freq[j]), fap


def ripiega(t, y, P, t0=0.0, nbin=100):
    """Curva ripiegata in nbin intervalli di fase (mediana per intervallo)."""
    fase = ((t - t0) / P) % 1.0
    idx = np.minimum((fase * nbin).astype(int), nbin - 1)
    med = np.full(nbin, np.nan)
    for k in range(nbin):
        v = y[idx == k]
        if len(v) >= 3:
            med[k] = np.median(v)
    if np.isnan(med).all():
        return med
    # buchi: interpolazione circolare
    x = np.arange(nbin)
    ok = ~np.isnan(med)
    med = np.interp(x, np.concatenate([x[ok] - nbin, x[ok], x[ok] + nbin]), np.tile(med[ok], 3))
    return med


# =============================================================================
# BRILLAMENTI
# =============================================================================
def brillamenti(t, f, cadenza, modello=None, con_indici=False):
    """Lampi brevi sopra la curva: [(tempo, ampiezza relativa, durata in minuti)].
    modello: la variazione periodica già trovata, da togliere prima di cercare i lampi."""
    if cadenza > 0.0105:            # più di 15 minuti: i brillamenti non si vedono bene
        return []
    if modello is not None:
        f = f - modello + 1.0
    fin = max(5, int(round(0.25 / cadenza)) | 1)   # mediana mobile di 6 ore
    from numpy.lib.stride_tricks import sliding_window_view
    pad = np.pad(f, fin // 2, mode="edge")
    tendenza = np.median(sliding_window_view(pad, fin), axis=1)[:len(f)]
    r = f - tendenza
    s = 1.4826 * np.median(np.abs(r - np.median(r)))
    if s <= 0:
        return []
    sopra = r > 3 * s
    out = []
    i, n = 0, len(r)
    min_punti = 3 if cadenza < 0.0035 else 2
    while i < n:
        if sopra[i]:
            j = i
            while j + 1 < n and sopra[j + 1] and t[j + 1] - t[j] < 3 * cadenza:
                j += 1
            pezzo = r[i:j + 1]
            if j - i + 1 >= min_punti and pezzo.max() > 5 * s and np.argmax(pezzo) <= (j - i) / 2 + 0.5:
                x = (float(t[i + int(np.argmax(pezzo))]), float(pezzo.max()), float((t[j] - t[i] + cadenza) * 1440))
                out.append(x + (i, j) if con_indici else x)
            i = j + 1
        else:
            i += 1
    return out


# =============================================================================
# CLASSIFICAZIONE
# =============================================================================
def mag(fl):
    return -2.5 * math.log10(max(fl, 1e-6))


def analizza(t, f, e, info, missione):
    """Restituisce None (niente di significativo) oppure un dict con periodo, tipo proposto, ecc."""
    cad = info.get("cadenza", 0.02)
    base = t.max() - t.min()
    sig = rumore_punti(f)
    # i brillamenti non devono disturbare il periodogramma: si tolgono solo i loro punti
    via = np.zeros(len(t), bool)
    for x in brillamenti(t, f, cad, con_indici=True):
        via[max(0, x[3] - 2):x[4] + 6] = True
    tc, fc = t[~via], f[~via]
    flare = []
    pmin = max(0.02, 4 * cad)
    pmax = min(base / 2.5, 45.0)
    risultato = None
    if len(tc) > 200 and pmax > pmin * 2:
        freq, pot, P, fap = periodogramma(tc, fc, pmin, pmax)
        b1 = ripiega(tc, fc, P)
        amp1 = float(np.nanmax(b1) - np.nanmin(b1))
        punti_bin = max(1, len(tc) / 100)
        snr = amp1 / (sig / math.sqrt(punti_bin)) if sig > 0 else 0
        sistematico = any(abs(P / q - 1) < 0.03 for q in PERIODI_SISTEMATICI) or P > 0.93 * pmax
        amp_mag = abs(mag(1 - amp1 / 2) - mag(1 + amp1 / 2))
        if fap < FAP_MAX and amp_mag >= AMPIEZZA_MIN_MAG and snr >= SNR_MIN and not sistematico:
            risultato = classifica(tc, fc, P, info, sig)
            risultato.update(fap=fap, snr=snr, freq=freq, pot=pot, sistematico=False)
    # brillamenti cercati dopo aver tolto la variazione periodica (le pulsanti salgono in fretta,
    # ma non sono brillamenti)
    if risultato is not None and risultato.get("periodo"):
        Pm = risultato["periodo"]
        b = ripiega(tc, fc, Pm, 0.0, 200)
        idx = np.minimum((((t / Pm) % 1.0) * 200).astype(int), 199)
        flare = brillamenti(t, f, cad, modello=b[idx])
        # i "lampi" che cadono dove la curva periodica cambia in fretta (salita di una RR Lyrae,
        # entrata in eclisse) sono resti del modello, non brillamenti veri
        g = np.abs(np.gradient(np.concatenate([b[-2:], b, b[:2]])))[2:-2]
        soglia = np.percentile(g, 75)
        flare = [x for x in flare if g[min(199, int(((x[0] / Pm) % 1.0) * 200))] <= soglia]
        if risultato["ampiezza"] and flare:
            flare = [x for x in flare if x[1] > 0.15 * risultato["ampiezza"]]
    else:
        flare = brillamenti(t, f, cad)
    if len(flare) >= 2:
        teff = info.get("teff")
        if risultato is None:
            risultato = dict(tipo="UV", periodo=None, epoca=None, ampiezza=None, durata=None, note=[],
                             max_rel=None, min_rel=None, fap=None, snr=None, freq=None, pot=None, sicuro=True)
        elif risultato["tipo"] in ("ROT", "BY"):
            risultato["tipo"] = risultato["tipo"] + "+UV"
        else:
            risultato["note"].append(f"{len(flare)} flares")
        if teff and teff > 5500 and risultato["tipo"] == "UV":
            risultato["note"].append("flares on a star hotter than 5500 K: check for contamination")
            risultato["sicuro"] = False
    if risultato is None:
        return None
    risultato["flare"] = flare
    risultato["rumore"] = sig
    risultato["base_giorni"] = base
    # magnitudini massima e minima nella banda della missione
    m0 = info.get("mag")
    if flare and risultato["tipo"].endswith("UV"):
        risultato["max_rel"] = max(risultato.get("max_rel") or 1.0, 1 + max(x[1] for x in flare))
        if risultato.get("min_rel") is None:
            risultato["min_rel"] = float(np.percentile(f, 1))
    if m0 is not None and risultato.get("max_rel"):
        risultato["mag_max"] = m0 + mag(risultato["max_rel"])
        risultato["mag_min"] = m0 + mag(risultato["min_rel"])
    crowd = info.get("crowdsap")
    if crowd is not None and crowd < CROWD_MIN:
        risultato["note"].append(f"CROWDSAP = {crowd:.2f}: the aperture is contaminated by nearby stars; "
                                 "the variability may belong to a neighbour")
        risultato["sicuro"] = False
    return risultato


def dispersione(t, f, P, nbin=60):
    """Quanto i punti si discostano dalla curva ripiegata media (più piccolo = periodo migliore)."""
    b = ripiega(t, f, P, 0.0, nbin)
    idx = np.minimum((((t / P) % 1.0) * nbin).astype(int), nbin - 1)
    return float(np.sqrt(np.nanmean((f - b[idx]) ** 2)))


def classifica(t, f, P, info, sig):
    """Tipo VSX proposto, periodo vero (P o 2P), epoca e durate, da regole sulla forma."""
    teff = info.get("teff") or float("nan")
    note = []
    # binarie a eclisse: il periodogramma spesso dà metà del periodo vero
    P2 = 2 * P
    b2 = ripiega(t, f, P2)
    n = len(b2)
    i1 = int(np.nanargmin(b2))
    i2 = (i1 + n // 2) % n
    vicino2 = [(i2 + k) % n for k in range(-6, 7)]
    d1 = 1 - b2[i1]
    d2 = 1 - np.nanmin(b2[vicino2])
    amp2 = float(np.nanmax(b2) - np.nanmin(b2))
    err_bin = sig / math.sqrt(max(1, len(t) / n))
    minimi_diversi = abs(d1 - d2) > max(5 * err_bin, 0.1 * amp2)
    # quanto della fase sta "in eclisse" (sotto metà ampiezza)
    livello = np.nanmax(b2) - amp2 / 2
    frazione_bassa = float(np.mean(b2 < livello))
    fuori = b2[b2 > np.nanmax(b2) - 0.25 * amp2]
    piatto_fuori = len(fuori) > 0.4 * n
    base = t.max() - t.min()
    raddoppiabile = P2 <= base / 2
    P_ls = P
    # il periodogramma trova spesso la metà del periodo vero: se ripiegando sul doppio i punti
    # si allineano molto meglio, il periodo vero è il doppio
    raddoppiato = raddoppiabile and dispersione(t, f, P2) < 0.8 * dispersione(t, f, P)
    due_minimi = d2 > 0.6 * d1 and frazione_bassa > 0.3          # due minimi simili per ciclo (contatto)
    if raddoppiato and not minimi_diversi and not due_minimi:
        P, P2, raddoppiabile = P2, 2 * P2, False       # una pulsante (o simile) al periodo doppio
        b2 = ripiega(t, f, P2)
    b1 = ripiega(t, f, P)
    amp1 = float(np.nanmax(b1) - np.nanmin(b1))
    amp_mag1 = abs(mag(np.nanmin(b1)) - mag(np.nanmax(b1)))
    # asimmetria: salita rapida e discesa lenta (pulsanti)
    imin, imax = int(np.nanargmin(b1)), int(np.nanargmax(b1))
    salita = ((imax - imin) % len(b1)) / len(b1)
    # le pulsanti tipo RR Lyrae salgono in fretta e scendono piano: le eclissi sono simmetriche
    pulsante = 0.2 <= P <= 1.2 and amp_mag1 > 0.15 and (salita < 0.3 or salita > 0.7)
    eclisse = raddoppiabile and not pulsante and (minimi_diversi or (piatto_fuori and frazione_bassa < 0.3)) and amp2 > 0.003
    # variazione continua con due minimi uguali, periodo doppio sotto 1 giorno: W UMa
    contatto = (raddoppiabile and not pulsante and not eclisse and P2 <= 1.0 and amp_mag1 >= 0.05
                and not (math.isfinite(teff) and teff >= 6800))
    if contatto:
        periodo, curva, tipo = P2, b2, "EW"
        im = int(np.nanargmin(curva))
        epoca_fase, durata = im / n, None
        massimi = [np.nanmax(b2[:n // 2]), np.nanmax(b2[n // 2:])]
        sicuro = abs(massimi[0] - massimi[1]) > 3 * err_bin       # massimi diversi (effetto O'Connell): più sicuro
    elif eclisse:
        periodo = P2
        if piatto_fuori and frazione_bassa < 0.3:
            tipo = "EA"
        elif P2 < 1.0 and not minimi_diversi:
            tipo = "EW"
        else:
            tipo = "EB" if minimi_diversi else ("EW" if P2 < 1.0 else "ELL")
        curva = b2
        im = int(np.nanargmin(curva))
        epoca_fase = im / n
        durata = None
        if tipo == "EA":       # durata dell'eclisse principale: intervalli contigui sotto il livello fuori eclisse
            soglia_ecl = np.nanmax(curva) - 0.1 * (np.nanmax(curva) - np.nanmin(curva))
            k = 1
            while k < n // 2 and curva[(im + k) % n] < soglia_ecl:
                k += 1
            h = 1
            while h < n // 2 and curva[(im - h) % n] < soglia_ecl:
                h += 1
            durata = 100 * (k + h - 1) / n
        sicuro = minimi_diversi or tipo == "EA"
    else:
        periodo = P
        curva = b1
        sicuro = True
        if pulsante and amp_mag1 > 0.2 and salita < 0.35:
            tipo = "RRAB"
        elif P < 0.25 and teff > 6300 if math.isfinite(teff) else P < 0.25:
            tipo = "HADS" if amp_mag1 > 0.3 else "DSCT"
        elif 0.25 <= P <= 1.2 and amp_mag1 > 0.2 and salita < 0.35:
            tipo = "RRAB"
        elif 0.2 <= P <= 0.5 and 0.15 < amp_mag1 < 0.7:
            tipo = "RRC"
        elif 0.3 <= P <= 3.5 and math.isfinite(teff) and 6700 <= teff <= 7600 and amp_mag1 < 0.1:
            tipo = "GDOR"
        elif P >= 0.3 and amp_mag1 < 0.2:
            tipo = "BY" if math.isfinite(teff) and teff < 5300 else "ROT"
            sicuro = math.isfinite(teff)
        else:
            tipo, sicuro = "VAR", False
        im = int(np.nanargmax(curva))
        epoca_fase = im / len(curva)
        durata = 100 * salita if tipo in ("RRAB", "RRC", "HADS", "DSCT") else None
    if not sicuro:
        note.append("type uncertain")
    fase_min = int(np.nanargmin(curva)) / len(curva)
    return dict(tipo=tipo, periodo=periodo, epoca=t.min() + ((epoca_fase - (t.min() / periodo) % 1.0) % 1.0) * periodo,
                epoca_min=t.min() + ((fase_min - (t.min() / periodo) % 1.0) % 1.0) * periodo,
                ampiezza=float(np.nanmax(curva) - np.nanmin(curva)), durata=durata, note=note,
                max_rel=float(np.nanmax(curva)), min_rel=float(np.nanmin(curva)), sicuro=sicuro, P_periodogramma=P_ls)


# =============================================================================
# CATALOGO VSX (copia HEASARC, aggiornata ogni settimana)
# =============================================================================
def cerca_vsx(ra, dec, raggio=RAGGIO_VSX_ARCSEC):
    """Variabili VSX entro il raggio: [(nome, tipo, periodo, separazione")]. Errori di rete: eccezione."""
    r = raggio / 3600
    dove = (f"CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {ra:.6f}, {dec:.6f}, {r:.6f})) = 1")
    try:
        righe = interroga_heasarc("SELECT name, ra, dec, variability_type, period, variability_flag "
                                  "FROM aavsovsx WHERE " + dove)
    except ValueError as e:
        if "column" not in str(e).lower() and "colonna" not in str(e).lower():
            raise
        righe = interroga_heasarc("SELECT * FROM aavsovsx WHERE " + dove)   # nomi delle colonne cambiati
    out = []
    for x in righe:
        try:
            ra2, dec2 = float(x.get("ra")), float(x.get("dec"))
        except (TypeError, ValueError):
            continue
        sep = math.degrees(math.acos(max(-1, min(1, math.sin(math.radians(dec)) * math.sin(math.radians(dec2)) +
                                                   math.cos(math.radians(dec)) * math.cos(math.radians(dec2)) *
                                                   math.cos(math.radians(ra - ra2)))))) * 3600
        if str(x.get("variability_flag", "0")).strip() in ("2", "3"):
            continue                               # "costante" o doppione: non conta come variabile nota
        out.append((x.get("name", ""), x.get("variability_type", ""), x.get("period", ""), sep))
    return out


# Formati di risposta provati in ordine (1.2: HEASARC non accetta più "votable/td").
# Quello che funziona viene ricordato per le stelle successive.
FORMATI_HEASARC = ["csv", "votable", "votable/td"]
_formato_ok = [None]


def interroga_heasarc(query):
    """Query ADQL al TAP di HEASARC: lista di dizionari {colonna minuscola: testo}."""
    formati = ([_formato_ok[0]] if _formato_ok[0] else []) + [f for f in FORMATI_HEASARC if f != _formato_ok[0]]
    ultimo = None
    for fmt in formati:
        par = {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": fmt, "RESPONSEFORMAT": fmt, "QUERY": query}
        dati = M.scarica(URL_HEASARC + "?" + urllib.parse.urlencode(par), 60)
        try:
            righe = leggi_risposta_tap(dati)
        except ValueError as e:
            ultimo = e
            if "format" in str(e).lower() or "formato" in str(e).lower():
                continue                           # formato rifiutato: provo il prossimo
            raise
        _formato_ok[0] = fmt
        return righe
    raise ultimo or ValueError("nessun formato accettato da HEASARC")


def leggi_risposta_tap(dati):
    """Riconosce da solo VOTable (TABLEDATA o binaria) o CSV."""
    inizio = dati[:3000].lstrip()
    if b"<VOTABLE" in inizio or inizio.startswith(b"<?xml"):
        return leggi_votable(dati)
    if inizio.startswith(b"<"):
        raise ValueError("risposta di HEASARC non valida (HTML?)")
    import io
    testo = dati.decode("utf-8", "replace")
    if "unsupported format" in testo[:500].lower():
        raise ValueError("errore di HEASARC: " + testo[:200])
    lettore = csv.reader(io.StringIO(testo))
    try:
        campi = [c.strip().lower() for c in next(lettore)]
    except StopIteration:
        return []
    if not campi or not any(campi):
        raise ValueError("risposta di HEASARC non valida (CSV vuoto)")
    return [dict(zip(campi, [v.strip() for v in riga])) for riga in lettore if riga]


def leggi_votable(dati):
    import xml.etree.ElementTree as ET
    import io
    if b"<VOTABLE" not in dati[:3000]:
        raise ValueError("risposta di HEASARC non valida")
    campi, righe, binaria = [], [], False
    for ev, el in ET.iterparse(io.BytesIO(dati), events=("end",)):
        tag = el.tag.rsplit("}", 1)[-1]
        if tag == "FIELD":
            campi.append((el.get("name") or "").lower())
        elif tag == "TR":
            righe.append(dict(zip(campi, [(td.text or "").strip() for td in el if td.tag.rsplit("}", 1)[-1] == "TD"])))
            el.clear()
        elif tag in ("BINARY", "BINARY2", "FITS"):
            binaria = True
        elif tag == "INFO" and el.get("value", "").upper() == "ERROR":
            raise ValueError("errore di HEASARC: " + (el.text or "")[:200])
    if binaria and not righe:                      # VOTable binaria: la legge astropy
        from astropy.io.votable import parse_single_table
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            tab = parse_single_table(io.BytesIO(dati)).to_table()
        righe = []
        for r in tab:
            d = {}
            for c in tab.colnames:
                v = r[c]
                if hasattr(v, "mask") and v.mask:
                    v = ""
                elif isinstance(v, bytes):
                    v = v.decode("utf-8", "replace")
                d[c.lower()] = str(v).strip()
            righe.append(d)
    return righe


# =============================================================================
# SCHEDE (in inglese) E MODULO VSX
# =============================================================================
def nome_stella(ident, info):
    return ident.replace("  ", " ")


# TT - UTC (secondi): serve per passare da BJD(TDB) all'HJD(UTC) che chiede VSX
TT_MENO_UTC = 69.184          # dal 2017 (TESS); per Kepler/K2 (2009-2018) si usa 66.184-69.184


def bjd_in_hjd(bjd, ra, dec):
    """BJD(TDB) -> HJD(UTC). La differenza tra baricentro e centro del Sole è al massimo
    ~4 secondi; TT-UTC ~69 secondi. Risultato preciso a qualche secondo (VSX usa 4 decimali)."""
    tt_utc = TT_MENO_UTC if bjd > 2457754.5 else 66.184 if bjd > 2454832.5 else 65.184
    return bjd - tt_utc / 86400.0


def dati_tic(info, ident):
    """Dal catalogo TIC (MAST): nomi in altri cataloghi, coordinate Gaia J2000, magnitudini V e T,
    stelle vicine che possono mescolare la luce nei pixel di TESS (21")."""
    from astroquery.mast import Catalogs
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    out = {"vicini": []}
    m = re.search(r"TIC\s*(\d+)", ident)
    righe = None
    if m:
        righe = Catalogs.query_criteria(catalog="Tic", ID=int(m.group(1)))
    if righe is None or len(righe) == 0:
        c = SkyCoord(info["ra"], info["dec"], unit="deg")
        righe = Catalogs.query_region(c, radius=5 * u.arcsec, catalog="Tic")
        if len(righe):
            righe.sort("dstArcSec")
    if righe is None or len(righe) == 0:
        return out
    r = righe[0]

    def val(k):
        try:
            v = r[k]
            if hasattr(v, "mask") and v.mask:
                return None
            v = v.item() if hasattr(v, "item") else v
            return None if v in ("", None) or (isinstance(v, float) and not math.isfinite(v)) else v
        except (KeyError, IndexError):
            return None
    out.update(tic=val("ID"), ra=val("ra"), dec=val("dec"), Vmag=val("Vmag"), Tmag=val("Tmag"), teff=val("Teff"),
               hip=val("HIP"), tyc=val("TYC"), twomass=val("TWOMASS"), gaia=val("GAIA"), contratio=val("contratio"))
    # stelle vicine entro 2 pixel di TESS
    if out.get("ra") is not None:
        c = SkyCoord(out["ra"], out["dec"], unit="deg")
        v = Catalogs.query_region(c, radius=42 * u.arcsec, catalog="Tic")
        for x in v:
            try:
                if str(x["ID"]) == str(out["tic"]) or x["Tmag"] is None:
                    continue
                dmag = float(x["Tmag"]) - float(out["Tmag"] or 99)
                if dmag < 4:                                   # abbastanza luminosa da contare
                    out["vicini"].append((f"TIC {x['ID']}", float(x["Tmag"]), float(x["dstArcSec"])))
            except (TypeError, ValueError, KeyError):
                continue
        out["vicini"].sort(key=lambda z: z[1])
    return out


def nome_vsx(ident, tic):
    """Nome principale come chiede VSX: prima HIP, poi Tycho-2, poi 2MASS, infine TIC."""
    if tic.get("hip"):
        return f"HIP {tic['hip']}"
    if tic.get("tyc"):
        return f"TYC {tic['tyc']}"
    if tic.get("twomass"):
        return f"2MASS J{tic['twomass']}"
    if tic.get("tic"):
        return f"TIC {tic['tic']}"
    return ident


def cifre_periodo(P, base):
    """Decimali significativi del periodo: incertezza circa 0,01 x P^2 / durata dei dati."""
    sigma = 0.01 * P * P / max(base, 1.0)       # errore di fase ~1 % diviso per il numero di cicli
    return int(min(7, max(2, math.ceil(-math.log10(max(sigma, 1e-9))))))


def campi_vsx(ident, missione, lotto, info, v, tic=None):
    """I campi del modulo "New star" di VSX, pronti da copiare (regole del manuale VSX)."""
    tic = tic or {}
    off = OFFSET_BJD.get(missione, 0.0)
    fonte = {"tess2": f"TESS 2-min cadence (SPOC PDCSAP), sector {lotto.lstrip('s0') or lotto}",
             "tessffi": f"TESS full-frame images (TESS-SPOC PDCSAP), sector {lotto.lstrip('s0') or lotto}",
             "kepler": "Kepler long-cadence PDCSAP, all quarters", "k2": f"K2 campaign {lotto} PDCSAP"}.get(missione, missione)
    ra = tic.get("ra", info.get("ra"))
    dec = tic.get("dec", info.get("dec"))
    tipo = v["tipo"] + ("" if v.get("sicuro", True) else ":")
    # epoca: minimo per binarie, ellissoidali e stelle che ruotano; massimo per le pulsanti
    al_minimo = v["tipo"].split("+")[0] in ("EA", "EB", "EW", "ELL", "ROT", "BY")
    epoca = None
    if v.get("epoca") is not None and v.get("periodo"):
        e = v["epoca"]
        if v["tipo"].split("+")[0] in ("ROT", "BY"):      # per le stelle che ruotano VSX vuole il minimo
            e = v.get("epoca_min", e)
        epoca = bjd_in_hjd(e + off, ra, dec)
    note = [f"Discovered in {fonte} light curves (MAST)."]
    altri = [x for x in (f"TIC {tic['tic']}" if tic.get("tic") else "",
                         f"2MASS J{tic['twomass']}" if tic.get("twomass") else "",
                         f"Gaia DR2 {tic['gaia']}" if tic.get("gaia") else "") if x and x != nome_vsx(ident, tic)]
    if altri:
        note.append("Other IDs: " + ", ".join(altri) + ".")
    if v.get("periodo"):
        note.append(f"Period from a Lomb-Scargle periodogram (false-alarm probability {v['fap']:.0e})"
                    + (", doubled because of unequal minima or a better phase fit" if v["periodo"] > 1.5 * v.get("P_periodogramma", v["periodo"]) else "") + ".")
    if v.get("flare"):
        note.append(f"{len(v['flare'])} flares in {v['base_giorni']:.1f} days.")
    teff = tic.get("teff") or info.get("teff")
    if teff:
        note.append(f"Teff = {float(teff):.0f} K (TIC).")
    # mescolamento della luce (blending): obbligatorio dirlo con pixel grandi come quelli di TESS
    vic = tic.get("vicini", [])
    if vic:
        n0 = vic[0]
        note.append(f"Possible blending: {len(vic)} star(s) within 42\" (2 TESS pixels) brighter than T+4; brightest "
                    f"{n0[0]} (T = {n0[1]:.2f}, {n0[2]:.0f}\"). No blending correction applied"
                    + (f"; CROWDSAP = {info['crowdsap']:.2f}" if info.get("crowdsap") is not None else "") + ".")
    else:
        note.append("No star brighter than T+4 within 42\" (2 TESS pixels)"
                    + (f"; CROWDSAP = {info['crowdsap']:.2f}" if info.get("crowdsap") is not None else "") + ".")
    if tic.get("Vmag"):
        note.append(f"Range given in the {BANDA_VSX.get(missione, '')} band (catalogue magnitude plus measured flux range); "
                    f"V = {float(tic['Vmag']):.2f} (TIC).")
    note += [n[0].upper() + n[1:] + "." for n in v.get("note", []) if n and n != "type uncertain"]
    amp = abs((v.get("mag_min") or 0) - (v.get("mag_max") or 0))
    dec_mag = 3 if amp < 0.05 else 2
    P = v.get("periodo")
    durata = ""
    if v.get("durata") and v["tipo"] == "EA":
        durata = f"{v['durata']:.0f}"
    elif v.get("durata") and v["tipo"] in ("RRAB", "HADS"):
        durata = f"{v['durata']:.0f}"
    c = {
        "Name": nome_vsx(ident, tic),
        "RA (J2000)": deg_hms(ra), "Dec (J2000)": deg_dms(dec),
        "RA (deg)": f"{ra:.6f}" if ra is not None else "", "Dec (deg)": f"{dec:+.6f}" if dec is not None else "",
        "Position source": "Gaia DR2 via TIC 8" if tic.get("ra") is not None else "TESS/Kepler FITS header",
        "Type": tipo,
        "Max": f"{v['mag_max']:.{dec_mag}f}" if v.get("mag_max") is not None else "",
        "Min": f"{v['mag_min']:.{dec_mag}f}" if v.get("mag_min") is not None else "",
        "Passband": BANDA_VSX.get(missione, ""),
        "Period (d)": f"{P:.{cifre_periodo(P, v['base_giorni'])}f}" if P else "",
        "Epoch (HJD)": f"{epoca:.4f}" if epoca is not None else "",
        "Epoch is": ("minimum" if al_minimo else "maximum") if epoca is not None else "",
        "Rise / eclipse duration (% of period)": durata,
        "Spectral type": "",
        "Discoverer": AUTORE,
        "Remarks": " ".join(note),
        "Attachments": "curva_luce.png (light curve), fase.png (phase plot)",
    }
    return c


def deg_hms(ra):
    if ra is None:
        return ""
    h = ra / 15
    hh, mm = int(h), int((h % 1) * 60)
    ss = ((h * 60) % 1) * 60
    return f"{hh:02d} {mm:02d} {ss:05.2f}"


def deg_dms(dec):
    if dec is None:
        return ""
    s = "-" if dec < 0 else "+"
    d = abs(dec)
    dd, mm = int(d), int((d % 1) * 60)
    ss = ((d * 60) % 1) * 60
    return f"{s}{dd:02d} {mm:02d} {ss:04.1f}"


def allegati(missione, t, f, v, campi, cartella):
    """I due allegati per VSX: curva di luce (luminosità verso l'alto) e curva ripiegata."""
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    m0 = None
    try:
        m0 = float(campi["Max"]) - mag(v["max_rel"]) if campi.get("Max") and v.get("max_rel") else None
    except (TypeError, ValueError):
        m0 = None
    y = (m0 + -2.5 * np.log10(np.clip(f, 1e-6, None))) if m0 is not None else f
    etich = f"{campi.get('Passband', '')} magnitude" if m0 is not None else "Normalized flux"
    off = OFFSET_BJD.get(missione, 0)
    fig = Figure(figsize=(8, 4.5), dpi=110)
    FigureCanvasAgg(fig)
    ax = fig.add_subplot(111)
    ax.plot(t + off - 2450000, y, ".", ms=1.5, color="#203060")
    if m0 is not None:
        ax.invert_yaxis()                       # luminosità che cresce verso l'alto
    ax.set_xlabel("BJD - 2450000")
    ax.set_ylabel(etich)
    ax.set_title(f"{campi['Name']} — light curve", fontsize=10)
    fig.tight_layout()
    fig.savefig(Path(cartella) / "curva_luce.png")
    if v.get("periodo"):
        P = v["periodo"]
        e0 = v["epoca_min"] if campi.get("Epoch is") == "minimum" else v["epoca"]
        fase = ((t - e0) / P) % 1.0
        fig = Figure(figsize=(8, 4.5), dpi=110)
        FigureCanvasAgg(fig)
        ax = fig.add_subplot(111)
        ax.plot(fase, y, ".", ms=1.5, color="#203060")
        ax.plot(fase + 1, y, ".", ms=1.5, color="#203060")
        if m0 is not None:
            ax.invert_yaxis()
        ax.set_xlim(0, 2)
        ax.set_xlabel(f"Phase (P = {campi['Period (d)']} d, epoch {campi['Epoch (HJD)']} HJD)")
        ax.set_ylabel(etich)
        ax.set_title(f"{campi['Name']} — phase plot ({campi['Type']})", fontsize=10)
        fig.tight_layout()
        fig.savefig(Path(cartella) / "fase.png")


def grafici(ident, missione, t, f, v, info, campi, cartella):
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.backends.backend_pdf import PdfPages
    cartella = Path(cartella)
    cartella.mkdir(parents=True, exist_ok=True)
    fig = Figure(figsize=(11.7, 8.3), dpi=120)
    FigureCanvasAgg(fig)
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 1.2, 1.1], hspace=0.6, wspace=0.25)
    ax = fig.add_subplot(gs[0, :])
    ax.plot(t, f, ".", ms=1, color="0.35")
    for x in v.get("flare", []):
        ax.axvline(x[0], color="#e04040", alpha=0.35, lw=1)
    ax.set_title("Light curve (normalized flux)" + (" — red lines: flares" if v.get("flare") else ""), fontsize=9)
    ax.set_xlabel("Time (BJD - %d)" % OFFSET_BJD.get(missione, 0), fontsize=8)
    if v.get("periodo"):
        P = v["periodo"]
        ax = fig.add_subplot(gs[1, 0])
        fase = (((t - v["epoca"]) / P) + 0.25) % 1.0 - 0.25
        ax.plot(fase, f, ".", ms=1, color="0.5", alpha=0.6)
        ax.plot(fase + 1, f, ".", ms=1, color="0.5", alpha=0.6)
        b = ripiega(t, f, P, v["epoca"] - 0.25 * P, 100)
        xb = (np.arange(100) + 0.5) / 100 - 0.25
        ax.plot(np.concatenate([xb, xb + 1]), np.tile(b, 2), "-", color="#d03030", lw=1.5)
        ax.set_xlim(-0.25, 1.25)
        ax.set_title(f"Phase plot, P = {P:.6f} d", fontsize=9)
        ax.set_xlabel("Phase", fontsize=8)
        ax = fig.add_subplot(gs[1, 1])
        ax.plot(1 / v["freq"], v["pot"], lw=0.6, color="#204080")
        ax.axvline(v.get("P_periodogramma", P), color="#d03030", alpha=0.4, lw=3)
        ax.set_xscale("log")
        ax.set_title("Lomb-Scargle periodogram", fontsize=9)
        ax.set_xlabel("Period (days)", fontsize=8)
    ax = fig.add_subplot(gs[2, :])
    ax.axis("off")
    righe = [f"{campi['Name']}  ·  {campi['Type']} ({NOMI_TIPI.get(v['tipo'].split('+')[0], v['tipo'])})",
             f"RA {campi['RA (J2000)']}   Dec {campi['Dec (J2000)']}   (J2000)",
             f"Range: {campi['Max']} – {campi['Min']} {campi['Passband']}" if campi["Max"] else "Range: n/a",
             f"Period: {campi['Period (d)'] or 'n/a'} d    Epoch ({campi['Epoch is'] or '-'}): {campi['Epoch (HJD)'] or 'n/a'} HJD",
             "Not in AAVSO VSX within %d arcsec (checked %s)." % (RAGGIO_VSX_ARCSEC, datetime.now().strftime("%Y-%m-%d"))]
    import textwrap
    righe += textwrap.wrap(campi["Remarks"], 125)
    for i, r in enumerate(righe[:12]):
        ax.text(0, 1 - i * 0.085, r, fontsize=8.3, transform=ax.transAxes, va="top", weight="bold" if i == 0 else "normal")
    fig.suptitle(f"{NOME_PROGRAMMA}: new variable star candidate {campi['Name']}  —  {COPYRIGHT}", fontsize=11)
    fig.savefig(cartella / "grafici.png")
    allegati(missione, t, f, v, campi, cartella)
    with PdfPages(cartella / "scheda.pdf") as pdf:
        pdf.savefig(fig)
    testo = ["AAVSO VSX - New star submission (copy each field into the New Star form at https://vsx.aavso.org)",
             "=" * 78]
    for k, val in campi.items():
        testo.append(f"{k}: {val}")
    testo += ["", "Before submitting:",
              "  - look at curva_luce.png and fase.png: the variation must be clear and repeated;",
              "  - read the blending remark: with a bright neighbour the variability may belong to it;",
              "  - attach curva_luce.png and fase.png (two files per save, as VSX asks);",
              "  - VSX accepts at most ONE new-star submission per user per day;",
              "  - incomplete submissions are deleted after 60 days."]
    (cartella / "vsx_submission.txt").write_text("\n".join(testo) + "\n", encoding="utf-8")


# =============================================================================
# UNA STELLA
# =============================================================================
def prepara(cartella, log):
    """Prima del giro: prova che il catalogo VSX risponda (una stella nota: delta Cephei)."""
    r = cerca_vsx(337.29277, 58.41520, 60)
    log(f"Catalogo AAVSO VSX (HEASARC) raggiungibile: {'sì' if r else 'risposta vuota'}")
    return None


def lavora_stella(missione, lotto, ident, sorgente, cartella, noti, n_thread=1, motore=None, log=print):
    cartella = Path(cartella)
    dl = cartella / M.CARTELLA_DOWNLOAD / re.sub(r"\W+", "_", ident)
    dl.mkdir(parents=True, exist_ok=True)
    t_start = time.time()
    try:
        urls = M.file_kepler(sorgente) if missione == "kepler" else [sorgente]
        percorsi = []
        for u in urls:
            nome = re.sub(r"[^\w.\-]", "_", urllib.parse.unquote(u.split("/")[-1].split("uri=")[-1]))
            percorsi.append(M.scarica_file(u, dl / nome))
        t, f, e, info = leggi(percorsi)
        v = analizza(t, f, e, info, missione)
        esito, candidati, vsx = "costante", [], ""
        if v is not None:
            noti_vsx = cerca_vsx(info["ra"], info["dec"]) if info.get("ra") is not None else []
            if noti_vsx:
                esito = "già in VSX"
                vsx = f"{noti_vsx[0][0]} ({noti_vsx[0][1]}, {noti_vsx[0][3]:.0f}\")"
            else:
                esito = "NUOVA"
                try:
                    tic = dati_tic(info, ident)
                except Exception as e:
                    if M.errore_di_rete(e):
                        raise
                    tic = {}
                campi = campi_vsx(ident, missione, lotto, info, v, tic)
                cart = cartella / CARTELLA_CANDIDATI / re.sub(r"\W+", "_", ident)
                grafici(ident, missione, t, f, v, info, campi, cart)
                with _scrittura:
                    M.aggiungi_csv(cartella / FILE_VSX, dict(campi))
                candidati.append(dict(stella=ident, missione=missione, lotto=lotto, campi=campi, v=v,
                                      cartella=f"{CARTELLA_CANDIDATI}/{cart.name}", qualita=qualita(v, tic, info)))
        riga = {"data": datetime.now().strftime("%Y-%m-%d %H:%M"), "missione": missione, "lotto": lotto, "stella": ident,
                "magnitudine": f"{info.get('mag', float('nan')):.2f}", "punti": len(t),
                "giorni": f"{t.max() - t.min():.1f}", "esito": esito,
                "tipo": v["tipo"] if v else "", "periodo": f"{v['periodo']:.6f}" if v and v.get("periodo") else "",
                "brillamenti": len(v["flare"]) if v else 0, "vsx": vsx, "secondi": f"{time.time() - t_start:.0f}"}
        return riga, candidati
    finally:
        shutil.rmtree(dl, ignore_errors=True)


def qualita(v, tic, info):
    """Punteggio 0-100: quali inviare per prime (VSX accetta una segnalazione al giorno)."""
    q = 50.0
    q += min(25, 4 * math.log10(max(v.get("snr") or 1, 1)))
    if v.get("sicuro", True):
        q += 15
    if tic.get("vicini"):
        q -= 10 * min(3, len(tic["vicini"]))
    c = info.get("crowdsap")
    if c is not None:
        q -= max(0, (0.95 - c) * 60)
    if v["tipo"] in ("EA", "EB", "EW", "RRAB", "RRC", "HADS"):
        q += 10                                   # tipi facili da confermare a occhio
    return round(max(0, min(100, q)))


def riga_candidato(c):
    r = {"data": datetime.now().strftime("%Y-%m-%d %H:%M"), "missione": c["missione"], "lotto": c["lotto"],
         "cartella": c["cartella"], "qualita": c.get("qualita", "")}
    r.update(c["campi"])
    return r


# =============================================================================
# PAGINA DEI RISULTATI (si apre dal PC con un doppio clic)
# =============================================================================
def scrivi_pagina(cartella):
    cartella = Path(cartella)
    p = cartella / FILE_CANDIDATI
    righe = []
    if p.exists():
        with open(p, encoding="utf-8-sig") as fh:
            righe = list(csv.DictReader(fh, delimiter=";"))
    # prima le più solide: VSX accetta una sola segnalazione al giorno
    def q(r):
        try:
            return float(r.get("qualita") or 0)
        except ValueError:
            return 0.0
    righe.sort(key=q, reverse=True)
    stato = {}
    try:
        stato = json.loads((cartella / "servizio_attivo.json").read_text(encoding="utf-8"))
    except (OSError, ValueError):
        pass
    tipi = {}
    for r in righe:
        tipi[r.get("Type", "?").rstrip(":")] = tipi.get(r.get("Type", "?").rstrip(":"), 0) + 1
    e = html.escape
    schede = []
    for r in righe:
        c = r.get("cartella", "")
        schede.append(f"""<div class="s"><a href="{e(c)}/grafici.png"><img src="{e(c)}/grafici.png" loading="lazy"></a>
<div><b>{e(r.get('Name', ''))}</b> <span class="t">{e(r.get('Type', ''))}</span> <span class="q">qualità {e(str(r.get('qualita', '')))}</span><br>
P = {e(r.get('Period (d)', '') or '–')} d · {e(r.get('Max', ''))}–{e(r.get('Min', ''))} {e(r.get('Passband', ''))}
· epoca {e(r.get('Epoch (HJD)', '') or '–')} HJD<br>
<small>{e(r.get('RA (J2000)', ''))} {e(r.get('Dec (J2000)', ''))} · {e(r.get('data', ''))}</small><br>
<a href="{e(c)}/vsx_submission.txt">modulo VSX</a> · <a href="{e(c)}/curva_luce.png">allegato 1</a> ·
<a href="{e(c)}/fase.png">allegato 2</a> · <a href="{e(c)}/scheda.pdf">scheda PDF</a></div></div>""")
    pagina = f"""<!doctype html><html lang="it"><head><meta charset="utf-8"><title>Variable Star Searcher</title>
<meta name="viewport" content="width=device-width, initial-scale=1"><meta http-equiv="refresh" content="600">
<style>body{{font:17px/1.4 system-ui,sans-serif;background:#0b0f19;color:#e6ecf8;margin:0;padding:16px}}
h1{{margin:0 0 4px}}.n{{color:#9aa8c4}}.g{{display:grid;grid-template-columns:repeat(auto-fill,minmax(340px,1fr));gap:12px;margin-top:14px}}
.s{{background:#141c2e;border:1px solid #26324d;border-radius:10px;padding:8px}}.s img{{width:100%;border-radius:6px;background:#fff}}
.t{{background:#1f4f8a;border-radius:6px;padding:0 6px;font-size:.85em}}.q{{color:#ffd36b;font-size:.85em}}a{{color:#8cc4ff}}</style></head><body>
<h1>⭐ Variable Star Searcher</h1>
<div class="n">Variabili nuove (non presenti in AAVSO VSX): <b>{len(righe)}</b> · {', '.join(f'{k} {n}' for k, n in sorted(tipi.items()))}<br>
Stato: {e(str(stato.get('stato', '–')))} · {e(str(stato.get('missione', '')))} {e(str(stato.get('lotto', '')))}
{stato.get('fatte', '')}/{stato.get('totale', '')} · aggiornata {datetime.now():%Y-%m-%d %H:%M} · {COPYRIGHT}</div>
<p class="n">Per inviarle: accedi a <a href="https://vsx.aavso.org">vsx.aavso.org</a>, apri il modulo per una nuova
variabile ("New Star"), copia i campi da <i>modulo VSX</i> e allega <i>allegato 1</i> e <i>allegato 2</i>.
VSX accetta <b>una sola segnalazione al giorno</b>: le variabili sono in ordine di qualità, dalla più solida.
Tutti i campi insieme: {FILE_VSX}.</p>
<div class="g">{''.join(schede) or '<p>Ancora nessuna variabile nuova.</p>'}</div></body></html>"""
    tmp = cartella / (FILE_PAGINA + ".tmp")
    tmp.write_text(pagina, encoding="utf-8")
    tmp.replace(cartella / FILE_PAGINA)


