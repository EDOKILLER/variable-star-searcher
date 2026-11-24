#!/bin/bash
# Variable Star Searcher - avvio del servizio nel container Docker del NAS
# © 2026 Edoardo Battistel
umask 000
cd /dati/programma || exit 1
export HOME=/dati/programma/_home MPLBACKEND=Agg PYTHONUNBUFFERED=1
mkdir -p "$HOME"
VENV=/dati/programma/_python
if [ ! -f "$VENV/.pronto_1.0" ]; then
    echo "Primo avvio: installo le librerie (5-10 minuti, serve Internet)…"
    until python -m venv "$VENV" && \
          "$VENV/bin/pip" install --no-cache-dir --upgrade pip && \
          "$VENV/bin/pip" install --no-cache-dir -r requisiti.txt; do
        echo "Installazione non riuscita (Internet assente?): riprovo tra 10 minuti"
        sleep 600
    done
    touch "$VENV/.pronto_1.0"
fi
exec nice -n 19 ionice -c 3 "$VENV/bin/python" servizio_var.py /dati
