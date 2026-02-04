#!/bin/bash
# stop_topkat.sh - Detener TopKAT Explorer App

APP_DIR="/home/jupyter-luisraul/topkat-explorer-app"
PID_FILE="$APP_DIR/topkat.pid"
LOG_FILE="$APP_DIR/topkat.log"

echo "========================================"
echo "DETENIENDO TOPKAT EXPLORER APP"
echo "Fecha: $(date)"
echo "========================================"

if [ -f "$PID_FILE" ]; then
    PID=$(cat "$PID_FILE")
    echo "Deteniendo proceso PID: $PID"
    
    # Intentar detener normalmente
    kill $PID 2>/dev/null
    sleep 3
    
    # Verificar si se detuvo
    if kill -0 $PID 2>/dev/null; then
        echo "Forzando detención..."
        kill -9 $PID
        sleep 1
    fi
    
    # Limpiar PID file
    rm "$PID_FILE"
    echo "✅ App detenida"
    
    # También matar cualquier proceso R suelto
    pkill -f "shiny::runApp" 2>/dev/null
else
    echo "⚠️  No se encontró PID file"
    echo "Buscando procesos Shiny..."
    
    PIDS=$(ps aux | grep "shiny::runApp" | grep -v grep | awk '{print $2}')
    if [ -n "$PIDS" ]; then
        echo "Encontrados procesos: $PIDS"
        kill $PIDS 2>/dev/null
        sleep 2
        kill -9 $PIDS 2>/dev/null
        echo "✅ Procesos eliminados"
    else
        echo "✅ No hay procesos Shiny corriendo"
    fi
fi

# Limpiar archivos temporales si existen
rm -f "$APP_DIR/.RData" 2>/dev/null
rm -f "$APP_DIR/.Rhistory" 2>/dev/null

echo ""
echo "📋 Para iniciar de nuevo: ./start_topkat.sh"
