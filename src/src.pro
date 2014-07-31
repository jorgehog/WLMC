include(../defaults.pri)

TEMPLATE = lib

TARGET = ../lib/WLMC

QMAKE_LFLAGS_DEBUG += -g

SOURCES += window.cpp \
    system.cpp \


HEADERS += window.h \
    system.h \
    windowparams.h

QMAKE_PRE_LINK += $(MKDIR) $$PWD/../lib $$shadowed($$PWD)/../lib

!equals(PWD, $${OUT_PWD}) {
    QMAKE_POST_LINK += $(COPY_DIR) $$OUT_PWD/../lib $$PWD/../
}
