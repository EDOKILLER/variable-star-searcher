# -*- coding: utf-8 -*-
"""
EXOPLANET SEARCHER - MOTORE DI CALCOLO
======================================
© 2026 Edoardo Battistel

Ricerca automatica di transiti planetari nelle curve di luce di TESS, Kepler e K2:
  1. elenco delle stelle (un settore TESS alla volta, poi le immagini a campo
     intero di TESS, poi Kepler e K2), dagli archivi MAST;
  2. download della curva di luce, pulizia (variazioni lente della stella, punti
     cattivi), ricerca dei transiti con TLS (o BLS se TLS non è disponibile);
  3. controlli contro i falsi positivi (pari/dispari, eclisse secondaria, forma
     a V, profondità da binaria) e confronto con i pianeti e candidati già noti
     (NASA Exoplanet Archive: pianeti confermati, TOI, KOI, K2);
  4. salvataggio di riassunto, candidati, grafici e scheda PDF; cancellazione
     dei FITS scaricati; si passa alla stella successiva.
"""

import csv
import io
import json
import math
import os
import re
import shutil
import sqlite3
import sys
import time
import traceback
import urllib.parse
import urllib.request
import warnings
from datetime import datetime
from pathlib import Path

import numpy as np

NOME_PROGRAMMA = "Exoplanet Searcher"
VERSIONE = "1.1"
AUTORE = "Edoardo Battistel"
COPYRIGHT = f"© 2026 {AUTORE}"
INTESTAZIONI_HTTP = {"User-Agent": "ExoplanetSearcher/1.1 (amateur astronomy software)"}

# --- archivi
URL_TESS_SCRIPT = "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_{s}_lc.sh"
URL_KEPLER = "https://archive.stsci.edu/missions/kepler/lightcurves/"
URL_K2 = "https://archive.stsci.edu/missions/k2/lightcurves/"
URL_MAST_FILE = "https://mast.stsci.edu/api/v0.1/Download/file?uri="
URL_EXOARCHIVE = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync"

MISSIONI = [("tess2", "TESS, curve a 2 minuti (SPOC)"),
            ("tessffi", "TESS, immagini a campo intero (TESS-SPOC)"),
            ("kepler", "Kepler (4 anni)"),
            ("k2", "K2 (campagne)")]

# --- soglie per i candidati
SDE_MIN = 9.0               # TLS: SDE 9 corrisponde a una probabilità di falso allarme ~0,01 %
SNR_MIN = 7.0
TRANSITI_MIN = 2
PROFONDITA_MAX = 0.05       # oltre il 5 % quasi sempre binarie a eclisse
RAGGIO_MAX_TERRE = 25.0     # più di ~2,2 raggi di Giove: non è un pianeta
MAX_SEGNALI = 3             # quanti segnali cercare in ogni stella (sistemi multipli)
PERIODO_MAX = 100.0         # giorni: oltre, la ricerca su 4 anni di Kepler diventerebbe lentissima
OVERSAMPLING_TLS = 3        # finezza della griglia dei periodi di TLS (più alto = più lento e più preciso)
TOLL_POSIZIONE_ARCSEC = 30.0
TOLL_PERIODO = 0.01         # 1 %

CARTELLA_DOWNLOAD = "_download"
CARTELLA_CANDIDATI = "candidati"
CARTELLA_CATALOGHI = "_cataloghi_noti"
FILE_STATO = "exo_stato.sqlite"
FILE_RIASSUNTO = "riassunto_stelle.csv"
FILE_CANDIDATI = "candidati.csv"
FILE_REGISTRO = "registro.txt"
FILE_SERVIZIO = "servizio_attivo.json"     # "battito" della ricerca sul NAS (lo legge il programma sul PC)

# Profili di ricerca: "normale" (PC) e "approfondito" (NAS: più lento, più accurato)
PROFILI = {
    "normale": dict(PERIODO_MAX=100.0, OVERSAMPLING_TLS=3, MAX_SEGNALI=3),
    "approfondito": dict(PERIODO_MAX=250.0, OVERSAMPLING_TLS=5, MAX_SEGNALI=5),
}


def applica_profilo(nome):
    """Imposta le costanti di ricerca del profilo scelto."""
    for k, v in PROFILI.get(nome, PROFILI["normale"]).items():
        globals()[k] = v


# =============================================================================
# RETE
# =============================================================================
def scarica(url, timeout=120, dati=None):
    req = urllib.request.Request(url, data=dati, headers=INTESTAZIONI_HTTP)
    with urllib.request.urlopen(req, timeout=timeout) as f:
        return f.read()


def scarica_file(url, dest, timeout=180):
    dest = Path(dest)
    tmp = dest.with_suffix(dest.suffix + ".part")
    req = urllib.request.Request(url, headers=INTESTAZIONI_HTTP)
    with urllib.request.urlopen(req, timeout=timeout) as f, open(tmp, "wb") as out:
        shutil.copyfileobj(f, out, 1 << 20)
    os.replace(tmp, dest)
    return dest


class ReteAssente(ConnectionError):
    """Internet non raggiungibile: si aspetta e si riprova, senza segnare errori."""


def errore_di_rete(e):
    """True se l'eccezione è un problema di connessione (non un file difettoso)."""
    import http.client
    import socket
    import urllib.error
    if isinstance(e, urllib.error.HTTPError):
        return e.code in (408, 429, 500, 502, 503, 504)
    if isinstance(e, (urllib.error.URLError, socket.timeout, TimeoutError, ConnectionError,
                      http.client.IncompleteRead, http.client.RemoteDisconnected, socket.gaierror)):
        return True
    testo = f"{type(e).__name__} {e}".lower()
    return any(x in testo for x in ("timed out", "timeout", "connection", "name resolution",
                                    "temporary failure", "network is unreachable", "remote end closed"))


def internet_disponibile(timeout=15):
    """Prova a raggiungere MAST e NASA Exoplanet Archive."""
    for url in ("https://archive.stsci.edu/", "https://exoplanetarchive.ipac.caltech.edu/"):
        try:
            req = urllib.request.Request(url, headers=INTESTAZIONI_HTTP, method="HEAD")
            with urllib.request.urlopen(req, timeout=timeout):
                return True
        except Exception as e:
            if not errore_di_rete(e):
                return True           # il server ha risposto (anche con un errore): la rete c'è
    return False


def elenco_html(url):
    """Nomi dei link di un elenco di cartella dell'archivio (Apache)."""
    testo = scarica(url, 60).decode("utf-8", "replace")
    return [x for x in re.findall(r'href="([^"?/][^"]*)"', testo) if not x.startswith(("http", "/", "#"))]


# =============================================================================
# ELENCHI DELLE STELLE
# =============================================================================
def ultimo_settore_tess(log=print, da=125):
    """Settore TESS più recente con le curve a 2 minuti già pubblicate."""
    if not internet_disponibile():
        raise ReteAssente("internet non disponibile")
    for s in range(da, 0, -1):
        try:
            req = urllib.request.Request(URL_TESS_SCRIPT.format(s=s), headers=INTESTAZIONI_HTTP, method="HEAD")
            with urllib.request.urlopen(req, timeout=20) as f:
                if f.status == 200:
                    return s
        except Exception:
            continue
    raise RuntimeError("nessun settore TESS trovato")


def stelle_tess_2min(settore):
    """[(id, url)] di tutte le curve a 2 minuti di un settore, dallo script curl di MAST."""
    testo = scarica(URL_TESS_SCRIPT.format(s=settore), 180).decode("utf-8", "replace")
    out = []
    for riga in testo.splitlines():
        m = re.search(r"-s(\d{4})-(\d{16})-\d{4}-[a-z]_lc\.fits\s+(https?://\S+)", riga)
        if m:
            out.append((f"TIC {int(m.group(2))}", m.group(3)))
    if not out:
        raise RuntimeError(f"settore {settore}: elenco vuoto")
    return out


def stelle_tess_ffi(settore, mag_max, log=print):
    """[(id, url)] delle curve TESS-SPOC (immagini a campo intero) fino a mag_max (Tmag)."""
    from astroquery.mast import Observations, Catalogs
    obs = Observations.query_criteria(provenance_name="TESS-SPOC", sequence_number=settore,
                                      dataproduct_type="timeseries")
    if len(obs) == 0:
        raise RuntimeError(f"TESS-SPOC settore {settore}: nessuna curva")
    tic = np.array([int(re.sub(r"\D", "", str(t)) or 0) for t in obs["target_name"]])
    log(f"TESS-SPOC settore {settore}: {len(obs)} stelle, lettura delle magnitudini...")
    mag = {}
    for i in range(0, len(tic), 800):
        blocco = [int(x) for x in tic[i:i + 800] if x]
        try:
            t = Catalogs.query_criteria(catalog="Tic", ID=blocco)
            for r in t:
                mag[int(r["ID"])] = float(r["Tmag"])
        except Exception as e:
            log(f"  magnitudini non lette per un blocco ({type(e).__name__})")
    scelte = [i for i, t in enumerate(tic) if mag.get(int(t), 99) <= mag_max]
    out = []
    for i in range(0, len(scelte), 200):
        sel = obs[scelte[i:i + 200]]
        prod = Observations.get_product_list(sel)
        for p in prod:
            nome = str(p["productFilename"])
            if nome.endswith("_lc.fits"):
                m = re.search(r"phot_0*(\d+)-s", nome)
                if m:
                    out.append((f"TIC {int(m.group(1))}", URL_MAST_FILE + str(p["dataURI"])))
    out.sort(key=lambda x: mag.get(int(x[0].split()[1]), 99))
    return out


def prefissi_kepler():
    return sorted(x.strip("/") for x in elenco_html(URL_KEPLER) if re.fullmatch(r"\d{4}/?", x))


def stelle_kepler(prefisso):
    """[(id, url_cartella)] delle stelle Kepler in una cartella dell'archivio (es. '0119')."""
    base = f"{URL_KEPLER}{prefisso}/"
    return [(f"KIC {int(x.strip('/'))}", base + x.strip("/") + "/")
            for x in elenco_html(base) if re.fullmatch(r"\d{9}/?", x)]


def campagne_k2():
    return sorted((x.strip("/") for x in elenco_html(URL_K2) if re.fullmatch(r"c\d{1,2}/?", x)),
                  key=lambda c: int(c[1:]))


def stelle_k2(campagna, log=print):
    """[(id, url_file)] delle curve K2 di una campagna (struttura c/N9/N5/ktwo..._llc.fits)."""
    out = []
    base = f"{URL_K2}{campagna}/"
    for a in elenco_html(base):
        if not re.fullmatch(r"\d{9}/?", a):
            continue
        for b in elenco_html(base + a.strip("/") + "/"):
            if not re.fullmatch(r"\d{5}/?", b):
                continue
            cart = base + a.strip("/") + "/" + b.strip("/") + "/"
            for f in elenco_html(cart):
                m = re.match(r"ktwo(\d{9})-c\d+_llc\.fits$", f)
                if m:
                    out.append((f"EPIC {int(m.group(1))}", cart + f))
    return out


def file_kepler(url_cartella):
    """URL dei file a cadenza lunga (llc) di una stella Kepler (tutti i trimestri)."""
    return [url_cartella + f for f in elenco_html(url_cartella) if f.endswith("_llc.fits")]


# =============================================================================
# LETTURA E PULIZIA DELLE CURVE DI LUCE
# =============================================================================
def leggi_curva(percorsi):
    """Unisce uno o più FITS di curve di luce (TESS, Kepler, K2, TESS-SPOC):
    tempo (giorni), flusso normalizzato, errore e informazioni sulla stella."""
    from astropy.io import fits
    T, F, E = [], [], []
    info = {}
    for p in percorsi:
        with fits.open(p, memmap=False) as h:
            hp = h[0].header
            if not info:
                for k_out, chiavi in (("ra", ("RA_OBJ",)), ("dec", ("DEC_OBJ",)),
                                      ("mag", ("TESSMAG", "KEPMAG")), ("teff", ("TEFF",)),
                                      ("raggio_stella", ("RADIUS",)), ("logg", ("LOGG",))):
                    for k in chiavi:
                        v = hp.get(k)
                        if isinstance(v, (int, float)) and math.isfinite(float(v)):
                            info[k_out] = float(v)
                            break
                info["missione"] = str(hp.get("MISSION", hp.get("TELESCOP", ""))).strip()
            d = h[1].data
            col = {c.upper(): c for c in d.columns.names}
            t = np.array(d[col["TIME"]], dtype=float)
            fl_col = col.get("PDCSAP_FLUX") or col.get("SAP_FLUX") or col.get("KSPSAP_FLUX")
            f = np.array(d[fl_col], dtype=float)
            e = np.array(d[col[fl_col.upper() + "_ERR"]], dtype=float) if fl_col.upper() + "_ERR" in col \
                else np.full_like(f, np.nan)
            q = np.array(d[col["QUALITY"]], dtype=int) if "QUALITY" in col else np.zeros(len(t), int)
            ok = np.isfinite(t) & np.isfinite(f) & (q == 0) & (f > 0)
            if ok.sum() < 50:
                continue
            med = np.median(f[ok])
            T.append(t[ok])
            F.append(f[ok] / med)
            E.append(e[ok] / med)
    if not T:
        raise ValueError("curva di luce vuota (nessun punto valido)")
    t, f, e = np.concatenate(T), np.concatenate(F), np.concatenate(E)
    o = np.argsort(t)
    t, f, e = t[o], f[o], e[o]
    if not np.isfinite(e).any() or np.nanmedian(e) <= 0:
        e = np.full_like(f, np.std(np.diff(f)) / math.sqrt(2))
    e = np.where(np.isfinite(e) & (e > 0), e, np.nanmedian(e))
    return t, f, e, info


def _segmenti(t, salto=0.5):
    tagli = np.where(np.diff(t) > salto)[0] + 1
    return np.split(np.arange(len(t)), tagli)


def pulisci(t, f, e, finestra_giorni=0.75):
    """Toglie le variazioni lente della stella (mediana mobile su ~0,75 giorni,
    calcolata per ogni tratto continuo) e i punti anomali verso l'alto (brillamenti).
    I transiti (cali brevi) restano."""
    from scipy.ndimage import median_filter
    tendenza = np.empty_like(f)
    cadenza = np.median(np.diff(t)) if len(t) > 1 else 0.02
    for idx in _segmenti(t):
        n = max(3, int(round(finestra_giorni / max(cadenza, 1e-4))) | 1)
        if len(idx) < n:
            tendenza[idx] = np.median(f[idx])
        else:
            tendenza[idx] = median_filter(f[idx], size=n, mode="nearest")
    p = f / tendenza
    mad = 1.4826 * np.median(np.abs(p - np.median(p)))
    ok = p < 1 + 3 * max(mad, 1e-6)                     # via i brillamenti
    return t[ok], p[ok], e[ok], tendenza


# =============================================================================
# RICERCA DEI TRANSITI
# =============================================================================
def motore_ricerca():
    """'TLS' se il pacchetto transitleastsquares funziona, altrimenti 'BLS'."""
    try:
        from transitleastsquares import transitleastsquares   # noqa: F401
        return "TLS"
    except Exception:
        return "BLS"


def _bls_numpy(t, y, periodi, durate):
    """BLS semplice in numpy (riserva per le prove senza astropy)."""
    w = np.ones_like(y)
    ym = y - np.average(y, weights=w)
    potenze, migliori = [], []
    for P in periodi:
        fase = (t % P) / P
        nb = 200
        b = np.minimum((fase * nb).astype(int), nb - 1)
        s = np.bincount(b, ym, nb)
        n = np.bincount(b, None, nb)
        best = (0, 0, 0)
        for d in durate:
            k = max(1, int(round(d / P * nb)))
            s2 = np.convolve(np.concatenate([s, s[:k]]), np.ones(k), "valid")[:nb]
            n2 = np.convolve(np.concatenate([n, n[:k]]), np.ones(k), "valid")[:nb]
            with np.errstate(invalid="ignore", divide="ignore"):
                pw = np.where(n2 > 2, s2 ** 2 / (n2 * (1 - n2 / len(y))), 0)
            pw = np.where(s2 < 0, pw, 0)
            j = int(np.argmax(pw))
            if pw[j] > best[0]:
                best = (pw[j], (j + k / 2) / nb * P, d)
        potenze.append(best[0])
        migliori.append(best)
    return np.array(potenze), migliori


def cerca_transiti(t, y, dy, raggio_stella=None, n_thread=1, motore=None):
    """Cerca il segnale di transito più forte. Restituisce un dict con periodo,
    T0, durata, profondità, SDE, SNR, numero di transiti, periodogramma, ecc."""
    base = t.max() - t.min()
    pmin, pmax = 0.5, max(1.0, min(base / 2.0, PERIODO_MAX))
    motore = motore or motore_ricerca()
    if motore == "TLS":
        from transitleastsquares import transitleastsquares
        modello = transitleastsquares(t, y, dy)
        argomenti = dict(period_min=pmin, period_max=pmax, oversampling_factor=OVERSAMPLING_TLS, duration_grid_step=1.1,
                         use_threads=max(1, int(n_thread)), show_progress_bar=False)
        if raggio_stella and 0.1 < raggio_stella < 10:
            argomenti.update(R_star=raggio_stella, R_star_min=max(0.08, raggio_stella * 0.5),
                             R_star_max=raggio_stella * 2)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            r = modello.power(**argomenti)
        if not math.isfinite(float(r.period)):
            return None
        return dict(motore="TLS", periodo=float(r.period), errore_periodo=float(r.period_uncertainty),
                    t0=float(r.T0), durata=float(r.duration), profondita=float(1 - r.depth),
                    sde=float(r.SDE), snr=float(r.snr), transiti=int(r.distinct_transit_count),
                    pari_dispari_sigma=float(r.odd_even_mismatch) if math.isfinite(float(r.odd_even_mismatch)) else 0.0,
                    periodi=np.asarray(r.periods), potenza=np.asarray(r.power))
    # BLS (astropy, o numpy come ultima riserva)
    durate = np.array([0.04, 0.06, 0.08, 0.12, 0.16, 0.24, 0.33])
    durate = durate[durate < pmin * 0.9]
    try:
        from astropy.timeseries import BoxLeastSquares
        bls = BoxLeastSquares(t, y, dy)
        res = bls.autopower(durate, minimum_period=pmin, maximum_period=pmax, frequency_factor=1.0,
                            objective="snr")
        periodi, potenza = np.asarray(res.period), np.asarray(res.power)
        j = int(np.argmax(potenza))
        P, t0, D = float(res.period[j]), float(res.transit_time[j]), float(res.duration[j])
    except ImportError:
        f = np.linspace(1 / pmax, 1 / pmin, 6000)
        periodi = 1 / f
        potenza, migliori = _bls_numpy(t, y, periodi, durate)
        j = int(np.argmax(potenza))
        P, t0, D = float(periodi[j]), float(migliori[j][1]), float(migliori[j][2])
    dentro = np.abs(((t - t0 + 0.5 * P) % P) - 0.5 * P) < D / 2
    if dentro.sum() < 3:
        return None
    prof = float(np.median(y[~dentro]) - np.mean(y[dentro]))
    rumore = 1.4826 * np.median(np.abs(y[~dentro] - np.median(y[~dentro])))
    snr = prof / (rumore / math.sqrt(dentro.sum())) if rumore > 0 else 0.0
    sde = float((potenza[j] - np.mean(potenza)) / (np.std(potenza) or 1))
    numeri = np.unique(np.round((t[dentro] - t0) / P))
    pari = dentro & (np.round((t - t0) / P) % 2 == 0)
    dispari = dentro & (np.round((t - t0) / P) % 2 == 1)
    pd = 0.0
    if pari.sum() > 2 and dispari.sum() > 2:
        dp = np.median(y[~dentro]) - np.mean(y[pari])
        dd = np.median(y[~dentro]) - np.mean(y[dispari])
        ee = rumore * math.sqrt(1 / pari.sum() + 1 / dispari.sum())
        pd = abs(dp - dd) / ee if ee > 0 else 0.0
    return dict(motore="BLS", periodo=P, errore_periodo=float("nan"), t0=t0, durata=D, profondita=prof,
                sde=sde, snr=float(snr), transiti=int(len(numeri)), pari_dispari_sigma=float(pd),
                periodi=periodi, potenza=potenza)


def maschera_transiti(t, periodo, t0, durata, fattore=1.5):
    return np.abs(((t - t0 + 0.5 * periodo) % periodo) - 0.5 * periodo) < fattore * durata / 2


def controlli(t, y, s, info):
    """Controlli contro i falsi positivi su un segnale trovato."""
    P, t0, D, prof = s["periodo"], s["t0"], s["durata"], s["profondita"]
    fase = ((t - t0) / P + 0.5) % 1 - 0.5
    dentro = np.abs(fase) < D / P / 2
    fuori = ~maschera_transiti(t, P, t0, D, 3)
    rumore = 1.4826 * np.median(np.abs(y[fuori] - np.median(y[fuori]))) if fuori.sum() > 10 else np.std(y)
    livello = np.median(y[fuori]) if fuori.sum() > 10 else 1.0
    # eclisse secondaria a metà periodo
    sec = np.abs(np.abs(fase) - 0.5) < D / P / 2
    prof_sec, sig_sec = 0.0, 0.0
    if sec.sum() > 3:
        prof_sec = float(livello - np.mean(y[sec]))
        sig_sec = prof_sec / (rumore / math.sqrt(sec.sum())) if rumore > 0 else 0.0
    # forma: profondità nella metà centrale / profondità totale (scatola ~1, V > ~1,5)
    centro = np.abs(fase) < D / P / 4
    forma = float((livello - np.mean(y[centro])) / prof) if centro.sum() > 2 and prof > 0 else float("nan")
    rs = info.get("raggio_stella") or 1.0
    raggio = math.sqrt(max(prof, 0)) * rs * 109.08
    avvisi, avvisi_en = [], []

    def avv(it, en):
        avvisi.append(it)
        avvisi_en.append(en)
    if s["sde"] < SDE_MIN:
        avv(f"SDE {s['sde']:.1f} < {SDE_MIN:g}", f"SDE {s['sde']:.1f} < {SDE_MIN:g}")
    if s["snr"] < SNR_MIN:
        avv(f"SNR {s['snr']:.1f} < {SNR_MIN:g}", f"SNR {s['snr']:.1f} < {SNR_MIN:g}")
    if s["transiti"] < TRANSITI_MIN:
        avv(f"solo {s['transiti']} transito", f"only {s['transiti']} transit")
    if prof > PROFONDITA_MAX:
        avv(f"profondità {prof * 100:.1f}% (binaria?)", f"depth {prof * 100:.1f}% (eclipsing binary?)")
    if raggio > RAGGIO_MAX_TERRE:
        avv(f"raggio {raggio:.0f} R⊕ troppo grande", f"radius {raggio:.0f} R⊕ too large")
    if s["pari_dispari_sigma"] > 3:
        avv(f"transiti pari/dispari diversi ({s['pari_dispari_sigma']:.1f}σ): binaria?",
            f"odd/even depth mismatch ({s['pari_dispari_sigma']:.1f}σ): eclipsing binary?")
    if sig_sec > 3:
        avv(f"eclisse secondaria ({sig_sec:.1f}σ): binaria?",
            f"secondary eclipse ({sig_sec:.1f}σ): eclipsing binary?")
    if math.isfinite(forma) and forma > 1.6:
        avv("forma a V (binaria radente?)", "V-shaped transit (grazing binary?)")
    return dict(raggio_terre=raggio, profondita_secondaria=prof_sec, sigma_secondaria=sig_sec,
                forma=forma, rumore=rumore, avvisi=avvisi, avvisi_en=avvisi_en,
                superato=not any(a.startswith(("SDE", "SNR", "solo", "profondità", "raggio")) for a in avvisi),
                pulito=not avvisi)


def analizza_stella(t, f, e, info, n_thread=1, motore=None, log=None):
    """Pulizia + ricerca di più segnali (dopo ogni segnale i suoi transiti vengono
    mascherati). Restituisce (curva pulita, lista dei segnali con i controlli)."""
    finestra = 0.75 if (t.max() - t.min()) < 200 else 1.0
    tp, yp, ep, tendenza = pulisci(t, f, e, finestra)
    segnali = []
    tt, yy, ee = tp, yp, ep
    for k in range(MAX_SEGNALI):
        if len(tt) < 200:
            break
        s = cerca_transiti(tt, yy, ee, info.get("raggio_stella"), n_thread, motore)
        if s is None:
            break
        c = controlli(tt, yy, s, info)
        s.update(c)
        s["numero"] = k + 1
        segnali.append(s)
        if not c["superato"]:
            break                           # segnale debole: inutile cercarne altri
        m = maschera_transiti(tt, s["periodo"], s["t0"], s["durata"], 2.0)
        tt, yy, ee = tt[~m], yy[~m], ee[~m]
    return (tp, yp, ep), segnali


# =============================================================================
# PIANETI E CANDIDATI GIÀ NOTI (NASA Exoplanet Archive)
# =============================================================================
QUERY_NOTI = {
    "confermati": "select pl_name,hostname,tic_id,ra,dec,pl_orbper,disc_facility from pscomppars",
    "toi": "select toi,tid,tfopwg_disp,pl_orbper,ra,dec from toi",
    "koi": "select kepid,kepoi_name,kepler_name,koi_disposition,koi_period,ra,dec from cumulative",
    "k2": "select pl_name,epic_hostname,disposition,pl_orbper,ra,dec from k2pandc",
}


def aggiorna_cataloghi_noti(cartella, log=print, giorni=7):
    """Scarica (al massimo una volta alla settimana) le tabelle dei pianeti e candidati noti."""
    c = Path(cartella) / CARTELLA_CATALOGHI
    c.mkdir(parents=True, exist_ok=True)
    for nome, q in QUERY_NOTI.items():
        p = c / f"{nome}.csv"
        if p.exists() and time.time() - p.stat().st_mtime < giorni * 86400:
            continue
        try:
            dati = scarica(URL_EXOARCHIVE + "?" + urllib.parse.urlencode({"query": q, "format": "csv"}), 300)
            if len(dati) > 100:
                p.write_bytes(dati)
                n_righe = dati.count(b"\n") - 1
                log(f"Catalogo dei noti aggiornato: {nome} ({n_righe} righe)")
        except Exception as e:
            log(f"Catalogo {nome} non aggiornato ({type(e).__name__}): uso la copia salvata" if p.exists()
                else f"Catalogo {nome} non disponibile ({type(e).__name__})")


class Noti:
    """Confronto di un segnale con pianeti confermati, TOI, KOI e candidati K2."""

    def __init__(self, cartella):
        self.righe = []
        c = Path(cartella) / CARTELLA_CATALOGHI
        for nome in QUERY_NOTI:
            p = c / f"{nome}.csv"
            if not p.exists():
                continue
            with open(p, encoding="utf-8", errors="replace") as f:
                for r in csv.DictReader(f):
                    try:
                        ra, dec = float(r.get("ra") or "nan"), float(r.get("dec") or "nan")
                    except ValueError:
                        continue
                    per = r.get("pl_orbper") or r.get("koi_period") or ""
                    try:
                        per = float(per)
                    except ValueError:
                        per = float("nan")
                    nomeogg = (r.get("pl_name") or (f"TOI-{r['toi']}" if r.get("toi") else "")
                               or r.get("kepler_name") or r.get("kepoi_name") or "")
                    stato = (r.get("tfopwg_disp") or r.get("koi_disposition") or r.get("disposition")
                             or ("CONFERMATO" if nome == "confermati" else ""))
                    ident = []
                    if r.get("tic_id"):
                        ident.append("TIC " + re.sub(r"\D", "", r["tic_id"]))
                    if r.get("tid"):
                        ident.append("TIC " + re.sub(r"\D", "", r["tid"]))
                    if r.get("kepid"):
                        ident.append("KIC " + re.sub(r"\D", "", r["kepid"]))
                    if r.get("epic_hostname"):
                        ident.append("EPIC " + re.sub(r"\D", "", r["epic_hostname"]))
                    self.righe.append(dict(catalogo=nome, nome=nomeogg, stato=stato, ra=ra, dec=dec,
                                           periodo=per, ident=set(i for i in ident if i.split()[-1])))
        self.ra = np.array([r["ra"] for r in self.righe]) if self.righe else np.zeros(0)
        self.dec = np.array([r["dec"] for r in self.righe]) if self.righe else np.zeros(0)

    def __len__(self):
        return len(self.righe)

    def confronta(self, ident, ra, dec, periodo):
        """('noto' | 'falso positivo noto' | 'nuovo in sistema noto' | 'nuovo', descrizione)"""
        vicini = [r for r in self.righe if ident in r["ident"]]
        if ra is not None and dec is not None and len(self.ra):
            d = np.degrees(np.arccos(np.clip(np.sin(np.radians(dec)) * np.sin(np.radians(self.dec)) +
                                             np.cos(np.radians(dec)) * np.cos(np.radians(self.dec)) *
                                             np.cos(np.radians(self.ra - ra)), -1, 1))) * 3600
            for i in np.where(d < TOLL_POSIZIONE_ARCSEC)[0]:
                if self.righe[i] not in vicini:
                    vicini.append(self.righe[i])
        for r in vicini:
            if math.isfinite(r["periodo"]) and r["periodo"] > 0:
                for mult in (1, 2, 0.5, 3, 1 / 3):
                    if abs(periodo / (r["periodo"] * mult) - 1) < TOLL_PERIODO:
                        fp = "FALSE" in r["stato"].upper() or r["stato"].upper() in ("FP", "FA")
                        tipo = "falso positivo noto" if fp else "noto"
                        rel = "" if mult == 1 else f" (periodo ×{mult:g})"
                        return tipo, f"{r['nome']} [{r['catalogo']}, {r['stato']}]{rel}"
        if vicini:
            return "nuovo in sistema noto", "stella con " + ", ".join(sorted({r["nome"] for r in vicini if r["nome"]})[:4])
        return "nuovo", ""


# =============================================================================
# STATO (riprendere da dove ci si era fermati)
# =============================================================================
class Stato:
    def __init__(self, cartella):
        self.cartella = Path(cartella)
        self.cartella.mkdir(parents=True, exist_ok=True)
        self.con = sqlite3.connect(self.cartella / FILE_STATO, check_same_thread=False, timeout=30)
        self.con.execute("CREATE TABLE IF NOT EXISTS stelle (missione TEXT, lotto TEXT, id TEXT, esito TEXT, "
                         "quando TEXT, PRIMARY KEY (missione, lotto, id))")
        self.con.execute("CREATE TABLE IF NOT EXISTS lotti (missione TEXT, lotto TEXT, finito INTEGER, "
                         "PRIMARY KEY (missione, lotto))")
        self.con.commit()

    def fatte(self, missione, lotto):
        return {r[0] for r in self.con.execute("SELECT id FROM stelle WHERE missione=? AND lotto=?",
                                               (missione, lotto))}

    def segna(self, missione, lotto, ident, esito):
        self.con.execute("INSERT OR REPLACE INTO stelle VALUES (?,?,?,?,?)",
                         (missione, lotto, ident, esito, datetime.now().strftime("%Y-%m-%d %H:%M")))
        self.con.commit()

    def lotto_finito(self, missione, lotto):
        r = self.con.execute("SELECT finito FROM lotti WHERE missione=? AND lotto=?", (missione, lotto)).fetchone()
        return bool(r and r[0])

    def segna_lotto(self, missione, lotto, finito=True):
        self.con.execute("INSERT OR REPLACE INTO lotti VALUES (?,?,?)", (missione, lotto, int(finito)))
        self.con.commit()

    def conteggi(self):
        tot = self.con.execute("SELECT COUNT(*) FROM stelle").fetchone()[0]
        per = dict(self.con.execute("SELECT missione, COUNT(*) FROM stelle GROUP BY missione").fetchall())
        return tot, per


def aggiungi_csv(percorso, riga):
    nuovo = not Path(percorso).exists()
    with open(percorso, "a", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=list(riga.keys()), delimiter=";")
        if nuovo:
            w.writeheader()
        w.writerow(riga)


# =============================================================================
# GRAFICI E SCHEDA DEL CANDIDATO
# =============================================================================
MISSIONI_EN = {"tess2": "TESS 2-min (SPOC)", "tessffi": "TESS FFI (TESS-SPOC)", "kepler": "Kepler", "k2": "K2"}
VERDETTI_EN = {"noto": "known", "falso positivo noto": "known false positive",
               "nuovo in sistema noto": "NEW in a known system", "nuovo": "NEW", "?": "not checked"}


def traduci_dettaglio(d):
    return (d.replace("stella con ", "— host star of ").replace("periodo ×", "period ×")
            .replace("CONFERMATO", "CONFIRMED").replace("confermati", "confirmed"))


def grafici_candidato(ident, missione, t_orig, f_orig, pulita, s, info, noto, cartella_out):
    """Grafici del candidato (PNG) e scheda riassuntiva (PDF di una pagina)."""
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.backends.backend_pdf import PdfPages
    t, y, _ = pulita
    P, t0, D = s["periodo"], s["t0"], s["durata"]
    cartella_out = Path(cartella_out)
    cartella_out.mkdir(parents=True, exist_ok=True)
    fig = Figure(figsize=(11.7, 8.3), dpi=120)
    FigureCanvasAgg(fig)
    gs = fig.add_gridspec(4, 2, height_ratios=[1, 1, 1.2, 1.4], hspace=0.75, wspace=0.25)
    ax = fig.add_subplot(gs[0, :])
    ax.plot(t_orig, f_orig, ".", ms=1, color="0.4")
    ax.set_title("Original light curve (normalized)", fontsize=9)
    ax.set_xlabel("Time (days)", fontsize=8)
    ax = fig.add_subplot(gs[1, :])
    ax.plot(t, y, ".", ms=1, color="0.3")
    m = maschera_transiti(t, P, t0, D, 1.0)
    ax.plot(t[m], y[m], ".", ms=2, color="#d03030")
    ax.set_title("Detrended light curve: detected transits in red", fontsize=9)
    ax.set_xlabel("Time (days)", fontsize=8)
    ax = fig.add_subplot(gs[2, 0])
    ax.plot(s["periodi"], s["potenza"], lw=0.6, color="#204080")
    ax.axvline(P, color="#d03030", alpha=0.4, lw=3)
    ax.set_xscale("log")
    ax.set_title(f"{s['motore']} periodogram (SDE {s['sde']:.1f})", fontsize=9)
    ax.set_xlabel("Period (days)", fontsize=8)
    fase = (((t - t0) / P + 0.5) % 1 - 0.5) * P * 24
    vis = np.abs(fase) < max(3 * D * 24, 3)
    ax = fig.add_subplot(gs[2, 1])
    ax.plot(fase[vis], y[vis], ".", ms=2, color="0.4")
    if vis.sum() > 20:
        b = np.linspace(fase[vis].min(), fase[vis].max(), 60)
        idx = np.digitize(fase[vis], b)
        med = [np.median(y[vis][idx == i]) if (idx == i).sum() > 2 else np.nan for i in range(1, len(b))]
        ax.plot((b[1:] + b[:-1]) / 2, med, "-", color="#d03030", lw=1.5)
    ax.set_title("Phase-folded transit", fontsize=9)
    ax.set_xlabel("Hours from mid-transit", fontsize=8)
    ax = fig.add_subplot(gs[3, 0])
    numero = np.round((t - t0) / P)
    for par, col, et in ((0, "#2060c0", "even"), (1, "#d07020", "odd")):
        sel = vis & (numero % 2 == par)
        ax.plot(fase[sel], y[sel], ".", ms=2, color=col, label=et, alpha=0.7)
    ax.legend(fontsize=7)
    ax.set_title(f"Odd and even transits (difference {s['pari_dispari_sigma']:.1f}σ)", fontsize=9)
    ax.set_xlabel("Hours from mid-transit", fontsize=8)
    ax = fig.add_subplot(gs[3, 1])
    ax.axis("off")
    avvisi = s.get("avvisi_en", s["avvisi"])
    righe = [f"{ident}  ·  {MISSIONI_EN.get(missione, missione)}", f"Period: {P:.5f} days", f"Duration: {D * 24:.2f} hours",
             f"Depth: {s['profondita'] * 1e6:.0f} ppm", f"Estimated radius: {s['raggio_terre']:.2f} R⊕"
             + ("" if info.get("raggio_stella") else " (assuming a Sun-like star)"),
             f"SDE {s['sde']:.1f}   SNR {s['snr']:.1f}   transits {s['transiti']}",
             f"Secondary eclipse: {s['sigma_secondaria']:.1f}σ   shape {s['forma']:.2f}",
             f"Catalogue check: {VERDETTI_EN.get(noto[0], noto[0])} {traduci_dettaglio(noto[1])}"[:90],
             "Warnings: " + ("; ".join(avvisi) if avvisi else "none")]
    if info.get("mag") is not None:
        righe.insert(1, f"Magnitude {info['mag']:.2f}   Teff {info.get('teff', float('nan')):.0f} K   "
                        f"R★ {info.get('raggio_stella', float('nan')):.2f} R☉")
    if info.get("ra") is not None and info.get("dec") is not None:
        righe.insert(2, f"RA {info['ra']:.5f}°   Dec {info['dec']:+.5f}°   (T0 = {t0:.5f})")
    for i, r in enumerate(righe):
        ax.text(0, 1 - i * 0.1, r, fontsize=8.2, transform=ax.transAxes, va="top",
                weight="bold" if i == 0 else "normal")
    fig.suptitle(f"{NOME_PROGRAMMA}: candidate {ident} #{s['numero']}  —  {COPYRIGHT}", fontsize=11)
    png = cartella_out / "grafici.png"
    fig.savefig(png)
    with PdfPages(cartella_out / "scheda.pdf") as pdf:
        pdf.savefig(fig)
    return png


# =============================================================================
# UNA STELLA DALL'INIZIO ALLA FINE
# =============================================================================
def lavora_stella(missione, lotto, ident, sorgente, cartella, noti, n_thread=1, motore=None, log=print):
    """Scarica, analizza, salva, cancella. Restituisce (riga del riassunto, candidati)."""
    cartella = Path(cartella)
    dl = cartella / CARTELLA_DOWNLOAD / re.sub(r"\W+", "_", ident)
    dl.mkdir(parents=True, exist_ok=True)
    t_start = time.time()
    try:
        if missione == "kepler":
            urls = file_kepler(sorgente)
        else:
            urls = [sorgente]
        percorsi = []
        for u in urls:
            nome = re.sub(r"[^\w.\-]", "_", urllib.parse.unquote(u.split("/")[-1].split("uri=")[-1]))
            percorsi.append(scarica_file(u, dl / nome))
        t, f, e, info = leggi_curva(percorsi)
        pulita, segnali = analizza_stella(t, f, e, info, n_thread, motore)
        candidati = []
        for s in segnali:
            if not s["superato"]:
                continue
            noto = noti.confronta(ident, info.get("ra"), info.get("dec"), s["periodo"]) if noti else ("?", "")
            s["noto"] = noto
            if noto[0] in ("nuovo", "nuovo in sistema noto"):
                sicuro = re.sub(r"\W+", "_", ident)
                c_out = cartella / CARTELLA_CANDIDATI / f"{sicuro}_{s['numero']}"
                grafici_candidato(ident, missione, t, f, pulita, s, info, noto, c_out)
                s["cartella"] = str(c_out)
                s["cartella_rel"] = f"{CARTELLA_CANDIDATI}/{sicuro}_{s['numero']}"
            candidati.append(s)
        migliore = segnali[0] if segnali else None
        riga = {"data": datetime.now().strftime("%Y-%m-%d %H:%M"), "missione": missione, "lotto": lotto,
                "stella": ident, "ra": f"{info.get('ra', float('nan')):.6f}", "dec": f"{info.get('dec', float('nan')):.6f}",
                "magnitudine": f"{info.get('mag', float('nan')):.2f}", "punti": len(t),
                "giorni": f"{t.max() - t.min():.1f}", "rumore_ppm": f"{(migliore or {}).get('rumore', float('nan')) * 1e6:.0f}",
                "motore": (migliore or {}).get("motore", ""),
                "periodo_migliore": f"{migliore['periodo']:.5f}" if migliore else "",
                "sde_migliore": f"{migliore['sde']:.2f}" if migliore else "",
                "candidati": len(candidati), "secondi": f"{time.time() - t_start:.0f}"}
        for s in candidati:
            s.update(stella=ident, missione=missione, lotto=lotto, info=info)
        return riga, candidati
    finally:
        shutil.rmtree(dl, ignore_errors=True)            # i FITS scaricati non servono più


def riga_candidato(s):
    info = s.get("info", {})
    return {"data": datetime.now().strftime("%Y-%m-%d %H:%M"), "stella": s["stella"], "segnale": s["numero"],
            "missione": s["missione"], "lotto": s["lotto"],
            "ra": f"{info.get('ra', float('nan')):.6f}", "dec": f"{info.get('dec', float('nan')):.6f}",
            "magnitudine": f"{info.get('mag', float('nan')):.2f}",
            "periodo_giorni": f"{s['periodo']:.5f}", "t0": f"{s['t0']:.5f}", "durata_ore": f"{s['durata'] * 24:.2f}",
            "profondita_ppm": f"{s['profondita'] * 1e6:.0f}", "raggio_terre": f"{s['raggio_terre']:.2f}",
            "sde": f"{s['sde']:.2f}", "snr": f"{s['snr']:.2f}", "transiti": s["transiti"],
            "pari_dispari_sigma": f"{s['pari_dispari_sigma']:.2f}", "secondaria_sigma": f"{s['sigma_secondaria']:.2f}",
            "forma": f"{s['forma']:.2f}", "controlli": "superati" if s["pulito"] else "; ".join(s["avvisi"]),
            "confronto": s["noto"][0], "dettaglio": s["noto"][1], "motore": s["motore"],
            "cartella": s.get("cartella_rel", s.get("cartella", ""))}


# =============================================================================
# LAVORO AUTOMATICO: missioni -> lotti -> stelle (usato dalla finestra e dal servizio NAS)
# =============================================================================
def n_nuclei():
    try:
        return max(1, len(os.sched_getaffinity(0)))
    except (AttributeError, OSError):
        return max(1, os.cpu_count() or 1)


class _CodaMuta:
    def put(self, x):
        pass


class Cacciatore:
    """Percorre missioni, lotti (settori, cartelle, campagne) e stelle.
    Se manca internet aspetta e riprova: le stelle non vengono segnate come errori."""

    ATTESA_MIN, ATTESA_MAX = 60, 1800          # secondi tra un tentativo e l'altro senza rete

    def __init__(self, cartella, opzioni, coda=None, ferma=None, thread_tls=None, lavoro=None, riga=None,
                 prepara=None, file_riassunto=None, file_candidati=None):
        """lavoro/riga/prepara/file_*: permettono di riusare il ciclo per altre ricerche
        (per esempio le stelle variabili) al posto di quella dei pianeti."""
        import threading
        self.lavoro = lavoro or lavora_stella
        self.riga = riga or riga_candidato
        self.prepara = prepara
        self.file_riassunto = file_riassunto or FILE_RIASSUNTO
        self.file_candidati = file_candidati or FILE_CANDIDATI
        self.cartella = Path(cartella)
        self.op = opzioni
        self.coda = coda or _CodaMuta()
        self.ferma = ferma or threading.Event()
        self.stato = Stato(self.cartella)
        self.motore = motore_ricerca()
        self.thread_tls = thread_tls
        self.senza_rete = False

    def log(self, m):
        self.coda.put(("log", m))
        try:
            with open(self.cartella / FILE_REGISTRO, "a", encoding="utf-8") as f:
                f.write(f"{datetime.now():%Y-%m-%d %H:%M:%S}  {m}\n")
        except OSError:
            pass

    def aspetta_rete(self, motivo=""):
        """Blocca finché internet non torna (o finché non viene chiesto di fermarsi)."""
        attesa = self.ATTESA_MIN
        self.senza_rete = True
        self.coda.put(("rete", False))
        self.log("Internet non disponibile" + (f" ({motivo})" if motivo else "") + ": aspetto e riprovo")
        while not self.ferma.is_set():
            if self.ferma.wait(attesa):
                break
            if internet_disponibile():
                self.log("Internet di nuovo disponibile: riprendo")
                break
            attesa = min(self.ATTESA_MAX, attesa * 2)
        self.senza_rete = False
        self.coda.put(("rete", True))

    def con_rete(self, funzione, descrizione):
        """Esegue funzione(); sui problemi di rete aspetta e riprova all'infinito."""
        while not self.ferma.is_set():
            try:
                return funzione()
            except Exception as e:
                if not errore_di_rete(e):
                    raise
                self.aspetta_rete(f"{descrizione}: {type(e).__name__}")
        return None

    def elenco(self, missione, lotto, funzione):
        """Elenco delle stelle di un lotto, salvato su disco per non riscaricarlo."""
        p = self.cartella / "_elenchi" / f"{missione}_{lotto}.json"
        if p.exists():
            return [tuple(x) for x in json.loads(p.read_text(encoding="utf-8"))]
        el = self.con_rete(funzione, f"elenco {missione} {lotto}")
        if el is None:
            return None
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps(el), encoding="utf-8")
        return el

    def lotti(self, missione):
        if missione in ("tess2", "tessffi"):
            ultimo = self.op.get("ultimo_settore") or self.con_rete(lambda: ultimo_settore_tess(self.log),
                                                                   "ultimo settore TESS")
            if not ultimo:
                return
            self.op["ultimo_settore"] = ultimo
            for s in range(ultimo, 0, -1):
                yield f"s{s:04d}", (lambda s=s: stelle_tess_2min(s)) if missione == "tess2" else \
                    (lambda s=s: stelle_tess_ffi(s, self.op.get("mag_ffi", 13.0), self.log))
        elif missione == "kepler":
            for pref in self.con_rete(prefissi_kepler, "cartelle Kepler") or []:
                yield pref, (lambda pref=pref: stelle_kepler(pref))
        elif missione == "k2":
            for c in self.con_rete(campagne_k2, "campagne K2") or []:
                yield c, (lambda c=c: stelle_k2(c, self.log))

    def esegui(self):
        self.log(f"Avvio: motore di ricerca {self.motore}"
                 + ("" if self.motore == "TLS" else " (TLS non disponibile: uso BLS)")
                 + f", periodo massimo {PERIODO_MAX:g} giorni, oversampling {OVERSAMPLING_TLS}, "
                   f"fino a {MAX_SEGNALI} segnali per stella")
        if self.prepara:
            noti = self.con_rete(self.prepara, "preparazione")
        else:
            self.con_rete(lambda: aggiorna_cataloghi_noti(self.cartella, self.log), "cataloghi dei noti")
            noti = Noti(self.cartella)
            self.log(f"Pianeti e candidati noti per il confronto: {len(noti)}")
        paralleli = max(1, int(self.op.get("paralleli", 2)))
        thread_tls = self.thread_tls or max(1, n_nuclei() // paralleli)
        for missione, etichetta in MISSIONI:
            if not self.op.get("missioni", {}).get(missione, True):
                continue
            try:
                for lotto, funzione in self.lotti(missione):
                    if self.ferma.is_set():
                        return
                    if self.stato.lotto_finito(missione, lotto):
                        continue
                    try:
                        stelle = self.elenco(missione, lotto, funzione)
                    except Exception as e:
                        self.log(f"{etichetta} {lotto}: elenco non disponibile ({type(e).__name__}: {str(e)[:120]})")
                        continue
                    if stelle is None:
                        return
                    fatte = self.stato.fatte(missione, lotto)
                    da_fare = [x for x in stelle if x[0] not in fatte]
                    self.log(f"{etichetta} {lotto}: {len(stelle)} stelle, da analizzare {len(da_fare)}")
                    self._lotto(missione, lotto, da_fare, len(stelle), noti, paralleli, thread_tls)
                    if self.ferma.is_set():
                        return
                    self.stato.segna_lotto(missione, lotto)
                    self.log(f"{etichetta} {lotto}: finito")
                    if not self.prepara:
                        try:
                            aggiorna_cataloghi_noti(self.cartella, self.log)
                        except Exception:
                            pass
            except Exception as e:
                if self.ferma.is_set():
                    return
                self.log(f"{etichetta}: ERRORE {type(e).__name__}: {e}")
        self.log("Tutte le missioni scelte sono state analizzate.")

    def _lotto(self, missione, lotto, da_fare, totale, noti, paralleli, thread_tls):
        from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
        fatte = totale - len(da_fare)
        t0, n0 = time.time(), fatte
        coda_stelle = list(da_fare)
        coda_stelle.reverse()                   # pop() dalla fine = ordine originale
        with ThreadPoolExecutor(max_workers=paralleli) as ex:
            in_corso = {}
            while True:
                while len(in_corso) < paralleli and coda_stelle and not self.ferma.is_set():
                    x = coda_stelle.pop()
                    fut = ex.submit(self.lavoro, missione, lotto, x[0], x[1], self.cartella, noti,
                                    thread_tls, self.motore, self.log)
                    in_corso[fut] = x
                if not in_corso:
                    break
                finiti, _ = wait(list(in_corso), timeout=5, return_when=FIRST_COMPLETED)
                rete_persa = False
                for fut in finiti:
                    x = in_corso.pop(fut)
                    ident = x[0]
                    try:
                        riga, candidati = fut.result()
                        aggiungi_csv(self.cartella / self.file_riassunto, riga)
                        for s in candidati:
                            r = self.riga(s)
                            aggiungi_csv(self.cartella / self.file_candidati, r)
                            self.coda.put(("candidato", r))
                        esito = f"{len(candidati)} candidati" if candidati else "nessun segnale"
                    except Exception as e:
                        if errore_di_rete(e):
                            coda_stelle.append(x)      # la riprendo quando torna la rete
                            rete_persa = True
                            continue
                        esito = f"errore: {type(e).__name__}: {str(e)[:100]}"
                        self.log(f"{ident}: {esito}")
                    self.stato.segna(missione, lotto, ident, esito)
                    fatte += 1
                    v = (fatte - n0) / max(1e-6, time.time() - t0)
                    resto = (totale - fatte) / v if v > 0 else 0
                    self.coda.put(("progresso", (missione, lotto, fatte, totale, v * 3600, resto, ident)))
                if rete_persa and not internet_disponibile():
                    # aspetto che finiscano le stelle in corso, poi aspetto la rete
                    wait(list(in_corso))
                    for fut in list(in_corso):
                        x = in_corso.pop(fut)
                        try:
                            riga, candidati = fut.result()
                        except Exception as e:
                            if errore_di_rete(e):
                                coda_stelle.append(x)
                                continue
                            self.stato.segna(missione, lotto, x[0], f"errore: {type(e).__name__}: {str(e)[:100]}")
                            continue
                        aggiungi_csv(self.cartella / self.file_riassunto, riga)
                        for s in candidati:
                            r = self.riga(s)
                            aggiungi_csv(self.cartella / self.file_candidati, r)
                            self.coda.put(("candidato", r))
                        self.stato.segna(missione, lotto, x[0],
                                         f"{len(candidati)} candidati" if candidati else "nessun segnale")
                        fatte += 1
                    self.aspetta_rete("download delle curve di luce")
                    t0, n0 = time.time(), fatte


# =============================================================================
# CURVA SINTETICA (per le prove)
# =============================================================================
def curva_sintetica(periodo=3.2, profondita=0.0015, durata=0.12, giorni=27.0, cadenza=2 / 1440, seme=1,
                    variabilita=0.003):
    rng = np.random.default_rng(seme)
    t = np.arange(0, giorni, cadenza) + 1500.0
    t = t[(t - t[0] < 13) | (t - t[0] > 14)]                 # buco a metà settore, come TESS
    f = 1 + variabilita * np.sin(2 * np.pi * (t - t[0]) / 5.3) + rng.normal(0, 0.0008, t.size)
    fase = ((t - 1500.7) / periodo + 0.5) % 1 - 0.5
    f[np.abs(fase) * periodo < durata / 2] -= profondita
    return t, f, np.full_like(f, 0.0008)
