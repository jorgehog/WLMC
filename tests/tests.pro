include(../defaults.pri)

TEMPLATE = app
CONFIG += console

INCLUDEPATH  += $$TOP_PWD/include

LIBS += -L../lib -lWLMC -lunittest++


TARGET = tests

SOURCES = testmain.cpp
