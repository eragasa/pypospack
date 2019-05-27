#!/usr/bin/env bash
# ADAPTED FROM https://wolfpaulus.com/technology/pythonlauncher/
start() {
  python `~/bin/pyposmat_slurm_daemon.py` &
  echo "pyposmat daemon started for python"
}

stop() {
      pid=`ps -ef | grep '[p]ython /etc/speech/TTSserver.py' | awk '{ print $2 }'`
      echo $pid
      kill $pid
      sleep 2
      echo "Server killed."
  }

  case "$1" in
    start)
      start
      ;;
    stop)
      stop
      ;;
    restart)
      stop
      start
      ;;
    *)
      echo "Usage: /etc/init.d/tornado-tts {start|stop|restart}"
      exit 1
  esac
  exit 0
}
