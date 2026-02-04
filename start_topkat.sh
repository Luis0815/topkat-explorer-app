#!/bin/bash
# start_topkat.sh - Iniciar TopKAT Explorer App

APP_DIR="/home/jupyter-luisraul/topkat-explorer-app"
PORT=3939
HOST="132.248.196.32"
LOG_FILE="$APP_DIR/topkat.log"
PID_FILE="$APP_DIR/topkat.pid"

echo "========================================"
echo "INICIANDO TOPKAT EXPLORER APP"
echo "Fecha: $(date)"
echo "IP: $HOST"
echo "Puerto: $PORT"
echo "========================================"

# Verificar si ya está corriendo
if [ -f "$PID_FILE" ]; then
    OLD_PID=$(cat "$PID_FILE")
    if kill -0 $OLD_PID 2>/dev/null; then
        echo "⚠️  La app YA ESTÁ CORRIENDO con PID: $OLD_PID"
        echo "   Para detener: ./stop_topkat.sh"
        echo "   Para ver logs: tail -f $LOG_FILE"
        exit 1
    else
        echo "Eliminando PID file obsoleto..."
        rm "$PID_FILE"
    fi
fi

# Cambiar al directorio
cd "$APP_DIR"

# Limpiar log anterior (opcional)
# echo "=== Nuevo inicio $(date) ===" > "$LOG_FILE"

# Iniciar app
echo "Iniciando Shiny app en $HOST:$PORT..."
echo "Esto puede tomar unos segundos..."

# Usando nohup para que sobreviva a desconexiones
nohup R -e "shiny::runApp(appDir = 'app.R', host = '$HOST', port = $PORT)" >> "$LOG_FILE" 2>&1 &
APP_PID=$!

# Guardar PID
echo $APP_PID > "$PID_FILE"
echo "✓ PID guardado: $APP_PID"

# Esperar que inicie
echo "Esperando que la app se inicie..."
sleep 8

# Verificar si está corriendo
if kill -0 $APP_PID 2>/dev/null; then
    echo "✅ APP INICIADA CORRECTAMENTE"
    echo ""
    echo "📊 INFORMACIÓN DE ACCESO:"
    echo "   URL: http://$HOST:$PORT"
    echo "   Puerto: $PORT"
    echo "   PID: $APP_PID"
    echo "   Logs: $LOG_FILE"
    echo ""
    echo "💡 COMANDOS ÚTILES:"
    echo "   Ver estado:   ./status_topkat.sh"
    echo "   Ver logs:     tail -f $LOG_FILE"
    echo "   Detener:      ./stop_topkat.sh"
    echo ""
    echo "🔒 La app seguirá corriendo aunque cierres esta terminal"
else
    echo "❌ ERROR: La app no se pudo iniciar"
    echo "Revisa los logs:"
    tail -20 "$LOG_FILE"
    rm "$PID_FILE"
    exit 1
fi
