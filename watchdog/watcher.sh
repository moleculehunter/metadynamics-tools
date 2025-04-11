#!/bin/bash

# === KONFIGURACJA ===
HILLS_FILE="hills/HILLS.8"
CHECK_INTERVAL=10                    # sekundy między sprawdzeniami
FREEZE_THRESHOLD=$((60))            # 10 min = 60 × 10s
MISSING_THRESHOLD=$((30))           # 5 min = 30 × 10s
SUBMIT_SCRIPT="submit.sh"
JOB_NAME="MDC40up"
SLURM_USER=$(whoami)
LOG_FILE="watch.log"

# === FUNKCJE ===

get_last_line() {
    tail -n 1 "$HILLS_FILE" 2>/dev/null
}

get_latest_jobid() {
    squeue -u "$SLURM_USER" --name="$JOB_NAME" --noheader --format="%i %V" \
    | sort -k2 | tail -n 1 | awk '{print $1}'
}

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG_FILE"
}

restart_job() {
    log "[RESTART $RESTART_COUNT] Killing frozen job: $JOBID"
    scancel "$JOBID"

    log "Waiting 60s before restart..."
    sleep 60

    log "Submitting new job..."
    sbatch "$SUBMIT_SCRIPT" >> "$LOG_FILE" 2>&1

    sleep 5
    NEW_JOBID=$(get_latest_jobid)

    if [[ -z "$NEW_JOBID" || ! "$NEW_JOBID" =~ ^[0-9]+$ ]]; then
        log "[ERROR] Could not retrieve new JOBID after restart!"
    else
        JOBID="$NEW_JOBID"
        log "New JOBID: $JOBID"
    fi
}

# === START ===

RESTART_COUNT=0
MISSING_COUNT=0

log "== Metadynamika Watcher Started =="

log "Waiting for initial job matching name: $JOB_NAME..."
while true; do
    JOBID=$(get_latest_jobid)
    if [[ -n "$JOBID" && "$JOBID" =~ ^[0-9]+$ ]]; then
        log "Initial JOBID: $JOBID"
        break
    fi
    sleep 5
done

# === GŁÓWNA PĘTLA ===
while true; do
    LAST_LINE=$(get_last_line)
    UNCHANGED_COUNT=0

    while true; do
        sleep "$CHECK_INTERVAL"

        # === WATCHDOG: brak jakiegokolwiek joba przez 5 min
        JOB_EXIST=$(squeue -u "$SLURM_USER" --name="$JOB_NAME" --noheader | wc -l)

        if (( JOB_EXIST == 0 )); then
            ((MISSING_COUNT++))
            echo "[DEBUG] No job found ($MISSING_COUNT/$MISSING_THRESHOLD)"
        else
            MISSING_COUNT=0
        fi

        if (( MISSING_COUNT >= MISSING_THRESHOLD )); then
            log "[WATCHDOG] Job '$JOB_NAME' not found for 5 minutes. Submitting new job."
            sbatch "$SUBMIT_SCRIPT" >> "$LOG_FILE" 2>&1
            sleep 5
            JOBID=$(get_latest_jobid)
            log "[WATCHDOG] New JOBID: $JOBID"
            MISSING_COUNT=0
            break
        fi

        # === MONITOR HILLS.8
        CURRENT_LINE=$(get_last_line)

        if [[ "$CURRENT_LINE" == "$LAST_LINE" ]]; then
            ((UNCHANGED_COUNT++))
            echo "[DEBUG] Line unchanged ($UNCHANGED_COUNT/$FREEZE_THRESHOLD)"
        else
            UNCHANGED_COUNT=0
            LAST_LINE="$CURRENT_LINE"
            echo "[DEBUG] Line changed. Reset counter."
        fi

        if (( UNCHANGED_COUNT >= FREEZE_THRESHOLD )); then
            JOB_STATE=$(squeue -j "$JOBID" --noheader --format="%T")

            if [[ "$JOB_STATE" == "RUNNING" ]]; then
                ((RESTART_COUNT++))
                log "[MONITOR] No updates to HILLS.8 for 10 min, job is RUNNING. Restarting."
                restart_job
                break
            else
                log "[MONITOR] No updates, but job state is '$JOB_STATE'. Skipping restart."
                UNCHANGED_COUNT=0
            fi
        fi
    done
done

