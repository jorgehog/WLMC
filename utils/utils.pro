TEMPLATE = subdirs
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += ordered

SUBDIRS += DCViz

HEADERS += lammpswriter/lammpswriter.h \
    BADAss/BADAss.h
