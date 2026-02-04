cat > monitor_topkat.sh << 'EOF'
#!/bin/bash
# monitor_topkat.sh - Monitorear y reiniciar si es necesario

APP_DIR="/home/jupyter-luisraul/topkat-explorer-app"
PID_FILE="$APP_DIR/topkat.pid"
LOG_FILE="$APP_DIR/monitor.log"
PORT=3939

echo "$(date): Iniciando verificación..." >> "$LOG_FILE"

# Verificar si el puerto está respondiendo
check_port() {
    timeout 5 bash -c "echo > /dev/tcp/132.248.196.32/$PORT" 2>/dev/null
    return $?
}

# Verificar proceso
if [ -f "$PID_FILE" ]; then
    PID=$(cat "$PID_FILE")
    if kill -0 $PID 2>/dev/null; then
        # Proceso existe, verificar puerto
        if check_port; then
            echo "$(date): OK - App corriendo y puerto activo" >> "$LOG_FILE"
            exit 0
        else
            echo "$(date): WARN - Proceso existe pero puerto no responde" >> "$LOG_FILE"
            # Reiniciar
            ./stop_topkat.sh >> "$LOG_FILE" 2>&1
            sleep 2
            ./start_topkat.sh >> "$LOG_FILE" 2>&1
        fi
    else
        echo "$(date): WARN - PID file existe pero proceso muerto" >> "$LOG_FILE"
        rm "$PID_FILE"
        ./start_topkat.sh >> "$LOG_FILE" 2>&1
    fi
else
    echo "$(date): WARN - No hay PID file, iniciando app" >> "$LOG_FILE"
    ./start_topkat.sh >> "$LOG_FILE" 2>&1
fi
EOF
