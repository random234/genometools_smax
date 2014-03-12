#!/bin/sh
clear
CC="ccache gcc" make curses=no with-sqlite=no cairo=no 64bit=no
