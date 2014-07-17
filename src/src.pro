include(../defaults.pri)

TEMPLATE = lib

TARGET = ../lib/WLMC

QMAKE_LFLAGS_DEBUG += -g

SOURCES += window.cpp \
    system.cpp \


HEADERS += window.h \
    system.h \
    windowparams.h
