#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
VARIABLE STAR SEARCHER - SERVIZIO PER IL NAS (senza finestra)
============================================================
© 2026 Edoardo Battistel

Gira di continuo dentro Docker sul NAS, a priorità minima, e cerca stelle variabili
non ancora catalogate (AAVSO VSX) nelle curve di luce di TESS, Kepler e K2:
  - legge le opzioni da var_config.json (nella cartella dei dati, creato al primo avvio);
  - se manca Internet aspetta e riprova (le stelle non vengono segnate come errori);
  - finito tutto, dorme 12 ore e ricontrolla se ci sono nuovi settori TESS;
  - scrive servizio_attivo.json ogni minuto e aggiorna risultati.html, la pagina
    con tutte le variabili nuove (si apre dal PC con un doppio clic).

Uso:  python servizio_var.py [cartella_dati]      (di base /dati)
      python servizio_var.py --prova               prova breve senza rete
"""

import json
import os
import signal
import sys
import threading
import time
import traceback
from datetime import datetime
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")
sys.path.insert(0, str(Path(__file__).resolve().parent))
import motore_exo as M  # noqa: E402
import motore_var as V  # noqa: E402

CONFIG_PREDEFINITA = {
    "_nota": "Opzioni di Variable Star Searcher sul NAS. Dopo una modifica riavvia il container.",
    "missioni": {"tess2": True, "tessffi": True, "kepler": True, "k2": True},
    "mag_ffi": 14.0,
    "paralleli": 3,
    "ore_di_pausa_a_fine_giro": 12,
}


class Battito:
    """Riceve gli eventi del Cacciatore e aggiorna servizio_attivo.json."""

    def __init__(self, cartella):
        self.p = Path(cartella) / M.FILE_SERVIZIO
        self.d = {"stato": "avvio", "rete": True, "fatte": 0, "totale": 0, "stelle_ora": 0.0,
                  "candidati_sessione": 0, "avvio": datetime.now().strftime("%Y-%m-%d %H:%M")}
        self.cartella = Path(cartella)
        self.lock = threading.Lock()
        self.pagina = True
        self.scrivi()

    def put(self, evento):
        tipo, d = evento
        with self.lock:
            if tipo == "log":
                print(f"{datetime.now():%Y-%m-%d %H:%M:%S}  {d}", flush=True)
                return
            if tipo == "progresso":
                mis, lotto, k, tot, ora, resto, ident = d
                self.d.update(stato="ricerca", missione=dict(M.MISSIONI).get(mis, mis), lotto=lotto,
                              fatte=k, totale=tot, stelle_ora=round(ora, 2), ore_alla_fine_lotto=round(resto / 3600, 1),
                              ultima=ident)
            elif tipo == "candidato":
                self.d["candidati_sessione"] += 1
                self.pagina = True
            elif tipo == "rete":
                self.d["rete"] = bool(d)
            elif tipo == "stato":
                self.d["stato"] = d
        self.scrivi()

    def scrivi(self):
        try:
            p_c = self.cartella / V.FILE_CANDIDATI
            n = 0
            if p_c.exists():
                with open(p_c, encoding="utf-8-sig") as f:
                    n = max(0, sum(1 for _ in f) - 1)
            self.d["candidati_totali"] = n
            self.d["ora"] = datetime.now().strftime("%Y-%m-%d %H:%M")
            tmp = self.p.with_suffix(".tmp")
            tmp.write_text(json.dumps(self.d, indent=1, ensure_ascii=False), encoding="utf-8")
            os.replace(tmp, self.p)
            if self.pagina or int(time.time()) % 600 < 60:
                self.pagina = False
                V.scrivi_pagina(self.cartella)
        except OSError:
            pass


def carica_config(cartella):
    p = Path(cartella) / "var_config.json"
    cfg = dict(CONFIG_PREDEFINITA)
    try:
        cfg.update(json.loads(p.read_text(encoding="utf-8")))
    except (OSError, ValueError):
        p.write_text(json.dumps(CONFIG_PREDEFINITA, indent=1, ensure_ascii=False), encoding="utf-8")
    return cfg


def abbassa_priorita():
    """Priorità minima anche se il container non usa nice/ionice."""
    try:
        os.nice(19 - os.nice(0))
    except OSError:
        pass
    try:
        import subprocess
        subprocess.run(["ionice", "-c", "3", "-p", str(os.getpid())], check=False,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except OSError:
        pass


def main(cartella):
    cartella = Path(cartella)
    cartella.mkdir(parents=True, exist_ok=True)
    abbassa_priorita()
    cfg = carica_config(cartella)
    ferma = threading.Event()

    def stop(*_):
        print("Arresto richiesto: finisco le stelle in corso…", flush=True)
        ferma.set()
    signal.signal(signal.SIGTERM, stop)
    signal.signal(signal.SIGINT, stop)

    battito = Battito(cartella)
    # battito anche quando non arrivano eventi (per esempio durante una stella Kepler lunga)
    def pulsa():
        while not ferma.wait(60):
            battito.scrivi()
    threading.Thread(target=pulsa, daemon=True).start()

    print(f"{V.NOME_PROGRAMMA} {V.VERSIONE} (servizio NAS) — {V.COPYRIGHT}", flush=True)
    print(f"Cartella: {cartella}  ·  nuclei visibili: {M.n_nuclei()}", flush=True)
    while not ferma.is_set():
        opzioni = {"missioni": cfg.get("missioni", {}), "mag_ffi": float(cfg.get("mag_ffi", 15.0)),
                   "paralleli": int(cfg.get("paralleli", 2))}
        cacc = M.Cacciatore(cartella, opzioni, battito, ferma, thread_tls=1, lavoro=V.lavora_stella,
                            riga=V.riga_candidato, prepara=lambda: V.prepara(cartella, cacc.log),
                            file_riassunto=V.FILE_RIASSUNTO, file_candidati=V.FILE_CANDIDATI)
        try:
            battito.put(("stato", "ricerca"))
            cacc.esegui()
        except Exception:
            cacc.log("ERRORE:\n" + traceback.format_exc())
            ferma.wait(600)              # errore imprevisto: riprovo tra 10 minuti
            continue
        if ferma.is_set():
            break
        ore = float(cfg.get("ore_di_pausa_a_fine_giro", 12))
        battito.put(("stato", f"in pausa fino a {datetime.fromtimestamp(time.time() + ore * 3600):%Y-%m-%d %H:%M}"))
        cacc.log(f"Giro completo: ricontrollo i nuovi settori tra {ore:g} ore")
        ferma.wait(ore * 3600)
    battito.put(("stato", "fermo"))
    try:
        battito.p.unlink()               # il PC smette subito di mostrare "in corso sul NAS"
    except OSError:
        pass
    print("Servizio fermato.", flush=True)


def prova():
    """Prova senza rete: classificazione su curve simulate, schede, modulo VSX, pagina, battito."""
    import tempfile
    import numpy as np
    d = Path(tempfile.mkdtemp(prefix="var_nas_prova_"))
    ok = True

    def v(c, testo):
        nonlocal ok
        print(("[ok] " if c else "[ERRORE] ") + testo)
        ok &= bool(c)
    rng = np.random.default_rng(3)
    cad = 10 / 1440
    t = np.arange(1500, 1527, cad)
    t = t[(t - 1500 < 13) | (t - 1500 > 14)]
    ph = ((t - 1500.3) / 2.7) % 1
    f = np.ones_like(t)
    f[np.abs(ph) < 0.03] -= 0.15
    f[np.abs(ph - 0.5) < 0.03] -= 0.05
    f += rng.normal(0, 0.0008, t.size)
    info = {"cadenza": cad, "mag": 10.0, "ra": 83.82, "dec": -5.39, "crowdsap": 0.95, "teff": 6000}
    r = V.analizza(t, f, np.full_like(f, 0.001), info, "tessffi")
    v(r is not None and r["tipo"] == "EA" and abs(r["periodo"] - 2.7) < 0.01, "binaria a eclisse simulata riconosciuta (EA, 2,7 giorni)")
    r2 = V.analizza(t, 1 + rng.normal(0, 0.0008, t.size), None, info, "tessffi")
    v(r2 is None, "stella costante: nessuna variabile")
    tic = {"tic": 123456789, "ra": 83.8221, "dec": -5.3911, "twomass": "05351730-0523280", "Vmag": 10.4, "Tmag": 10.0,
           "vicini": [("TIC 987", 11.2, 25.0)]}
    campi = V.campi_vsx("TIC 123456789", "tessffi", "s0050", info, r, tic)
    v(campi["Type"] == "EA" and campi["Epoch is"] == "minimum" and campi["Passband"] == "TESS"
      and campi["Name"].startswith("2MASS J") and "Possible blending" in campi["Remarks"]
      and campi["Epoch (HJD)"].count(".") == 1 and len(campi["Epoch (HJD)"].split(".")[1]) == 4, "campi del modulo VSX (regole del manuale)")
    V.grafici("TIC 123456789", "tessffi", t, f, r, info, campi, d / "variabili" / "TIC_123456789")
    cv = d / "variabili" / "TIC_123456789"
    v((cv / "vsx_submission.txt").exists() and (cv / "curva_luce.png").exists() and (cv / "fase.png").exists(),
      "scheda, allegati per VSX e modulo scritti")
    M.aggiungi_csv(d / V.FILE_CANDIDATI, V.riga_candidato(dict(stella="TIC 123456789", missione="tessffi", lotto="s0050",
                                                                campi=campi, cartella="variabili/TIC_123456789")))
    b = Battito(d)
    b.put(("progresso", ("tessffi", "s0050", 5, 100, 30.0, 3600, "TIC 1")))
    v((d / V.FILE_PAGINA).exists() and campi["Name"] in (d / V.FILE_PAGINA).read_text(encoding="utf-8"),
      "pagina risultati.html")
    import urllib.error
    v(M.errore_di_rete(urllib.error.URLError("no route")), "riconosce i problemi di rete")
    print(("PROVA SUPERATA: " if ok else "PROVA FALLITA: ") + str(d))
    return ok


if __name__ == "__main__":
    if "--prova" in sys.argv:
        sys.exit(0 if prova() else 1)
    else:
        argomenti = [a for a in sys.argv[1:] if not a.startswith("-")]
        main(argomenti[0] if argomenti else os.environ.get("EXO_DATI", "/dati"))
