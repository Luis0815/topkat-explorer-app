#!/bin/bash
# status_topkat.sh - Ver estado de TopKAT Explorer App

APP_DIR="/home/jupyter-luisraul/topkat-explorer-app"
PID_FILE="$APP_DIR/topkat.pid"
LOG_FILE="$APP_DIR/topkat.log"
HOST="132.248.196.32"
PORT=3939

echo "========================================"
echo "ESTADO DE TOPKAT EXPLORER APP"
echo "Fecha: $(date)"
echo "========================================"

# 1. Verificar por PID file
if [ -f "$PID_FILE" ]; then
    PID=$(cat "$PID_FILE")
    if kill -0 $PID 2>/dev/null; then
        echo "✅ APP CORRIENDO"
        echo "   PID: $PID"
        echo "   URL: http://$HOST:$PORT"
        echo "   Puerto: $PORT"
        
        # Obtener información del proceso
        echo ""
        echo "📊 INFORMACIÓN DEL PROCESO:"
        ps -p $PID -o pid,user,pcpu,pmem,etime,cmd | tail -1
        
        # Ver conexiones
        echo ""
        echo "🔗 CONEXIONES AL PUERTO $PORT:"
        netstat -tulpn 2>/dev/null | grep ":$PORT" || echo "   No se pudo verificar conexiones"
        
    else
        echo "❌ PID file existe pero proceso NO ACTIVO"
        echo "   PID: $PID"
        rm "$PID_FILE"
    fi
else
    echo "❌ NO HAY PID FILE"
fi

# 2. Verificar por nombre de proceso
echo ""
echo "🔍 BUSCANDO PROCESOS SHINY:"
SHINY_PROCS=$(ps aux | grep -E "shiny::runApp|runApp.*$PORT" | grep -v grep)
if [ -n "$SHINY_PROCS" ]; then
    echo "Procesos encontrados:"
    echo "$SHINY_PROCS"
else
    echo "   No hay procesos Shiny corriendo"
fi

# 3. Verificar si el puerto está en uso
echo ""
echo "📡 VERIFICANDO PUERTO $PORT:"
if command -v lsof >/dev/null 2>&1; then
    lsof -i :$PORT
elif command -v netstat >/dev/null 2>&1; then
    netstat -tulpn 2>/dev/null | grep ":$PORT"
else
    echo "   No se pudo verificar el puerto"
fi

# 4. Mostrar últimos logs
echo ""
echo "📝 ÚLTIMOS 10 REGISTROS DEL LOG:"
if [ -f "$LOG_FILE" ]; then
    tail -10 "$LOG_FILE"
else
    echo "   No hay archivo de log"
fi

echo ""
echo "========================================"
echo "💡 COMANDOS:"
echo "   Iniciar:  ./start_topkat.sh"
echo "   Detener:  ./stop_topkat.sh"
echo "   Ver logs: tail -f $LOG_FILE"
