CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CXX = gcc

!noccache {
    QMAKE_CXX = ccache $$QMAKE_CXX
}

COMMON_CXXFLAGS = -std=c++11

QMAKE_CXXFLAGS += \
    $$COMMON_CXXFLAGS

QMAKE_CXXFLAGS_DEBUG += \
    $$COMMON_CXXFLAGS

QMAKE_CXXFLAGS_RELEASE += \
    $$COMMON_CXXFLAGS \
    -O3 \
    -DNDEBUG \
    -DARMA_NO_DEBUG

QMAKE_LFLAGS_RELEASE -= -O1
QMAKE_LFLAGS_RELEASE += -O3

DEFINES += \
    ARMA_MAT_PREALLOC=3

INCLUDEPATH += $$PWD/utils

LIBS += -larmadillo -llapack -lblas


TOP_PWD = $$PWD
