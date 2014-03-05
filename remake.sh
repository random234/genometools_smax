#!/bin/sh
clear
CC="ccache gcc" make curses=no cairo=yes noopt=yes
