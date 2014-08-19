TEMPLATE = subdirs
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += ordered

SUBDIRS += src tests

OTHER_FILES += include/WLMC.h $(HOME)/tmp/WLC/wl-w.c
