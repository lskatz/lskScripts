#!/bin/sh

# Borrowed from 
# https://unix.stackexchange.com/questions/294299/how-to-renice-all-threads-and-children-of-one-process-on-linux

if [ "$#" -eq 0 ]; then
  echo "Finds all children pids of a process"
  echo "  Usage: $(basename $0) 012345"
  exit 1
fi

PID_LIST=
findpids() {
        for pid in /proc/$1/task/* ; do
                pid="$(basename "$pid")"
                PID_LIST="$PID_LIST$pid "
                if [ ! -e "/proc/$1/task/$pid/children" ]; then
                  continue;
                fi
                for cpid in $(cat /proc/$1/task/$pid/children 2>/dev/null) ; do
                        findpids $cpid
                done
        done
}

for pid in $@; do 

  if [ ! -e "/proc/$1/task" ]; then
    echo "ERROR: could not find pid $pid in the process list";
    exit 1
  fi

  findpids $1
done

echo $PID_LIST
