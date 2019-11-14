
#
#    use as
#
#    import RI.py
#
#    WLmin = RI.WLmin
#    WLmax = RI.WLmax
#    WLperiod = RI.WLperiod
#    WLx = RI.WLx
#    NumWLx = RI.NumWLx
#
#    RI.RIAul[i]   RI,RIAgl[i]   RI.epAu[i]   RI.epAg[i]
#

from scipy import array, interpolate, arange, zeros

RIAu = array([
    [292.4, 1.49, 1.878], [300.9, 1.53, 1.889], [310.7, 1.53, 1.893],
    [320.4, 1.54, 1.898], [331.5, 1.48, 1.883], [342.5, 1.48, 1.871],
    [354.2, 1.50, 1.866], [367.9, 1.48, 1.895], [381.5, 1.46, 1.933],
    [397.4, 1.47, 1.952], [413.3, 1.46, 1.958], [430.5, 1.45, 1.948],
    [450.9, 1.38, 1.914], [471.4, 1.31, 1.849], [495.9, 1.04, 1.833],
    [520.9, 0.62, 2.081], [548.6, 0.43, 2.455], [582.1, 0.29, 2.863],
    [616.8, 0.21, 3.272], [659.5, 0.14, 3.697], [704.5, 0.13, 4.103],
    [756.0, 0.14, 4.542], [821.1, 0.16, 5.083], [892.0, 0.17, 5.663],
    [984.0, 0.22, 6.350], [1088.0, 0.27, 7.150]])

RIAg = array([
    [292.4, 1.39, 1.161], [300.9, 1.34, 0.964], [310.7, 1.13, 0.616],
    [320.4, 0.81, 0.392], [331.5, 0.17, 0.829], [342.5, 0.14, 1.142],
    [354.2, 0.10, 1.419], [367.9, 0.07, 1.657], [381.5, 0.05, 1.864],
    [397.4, 0.05, 2.070], [413.3, 0.05, 2.275], [430.5, 0.04, 2.462],
    [450.9, 0.04, 2.657], [471.4, 0.05, 2.869], [495.9, 0.05, 3.093],
    [520.9, 0.05, 3.324], [548.6, 0.06, 3.586], [582.1, 0.05, 3.858],
    [616.8, 0.06, 4.152], [659.5, 0.05, 4.483], [704.5, 0.04, 4.838],
    [756.0, 0.03, 5.242], [821.1, 0.04, 5.727], [892.0, 0.04, 6.312],
    [984.0, 0.04, 6.992], [1088.0, 0.04, 7.795]])

NumWL = 26
WL = zeros(NumWL, dtype=int)
RIAuRe = zeros(NumWL, dtype=float)
RIAuIm = zeros(NumWL, dtype=float)
RIAgRe = zeros(NumWL, dtype=float)
RIAgIm = zeros(NumWL, dtype=float)

WLmin = 300
WLmax = 1000
WLperiod = 1
WLx = arange(WLmin, WLmax + 1, WLperiod)
NumWLx = int((WLmax + 1 - WLmin) / WLperiod)

for i in range(NumWL):
    WL[i] = RIAu[i, 0]
    RIAuRe[i] = RIAu[i, 1]
    RIAuIm[i] = RIAu[i, 2]
    RIAgRe[i] = RIAg[i, 1]
    RIAgIm[i] = RIAg[i, 2]

fRIAuReInt2 = interpolate.splrep(WL, RIAuRe, s=0)
RIAuReInt2 = interpolate.splev(WLx, fRIAuReInt2, der=0)

fRIAuImInt2 = interpolate.splrep(WL, RIAuIm, s=0)
RIAuImInt2 = interpolate.splev(WLx, fRIAuImInt2, der=0)

fRIAgReInt2 = interpolate.splrep(WL, RIAgRe, s=0)
RIAgReInt2 = interpolate.splev(WLx, fRIAgReInt2, der=0)

fRIAgImInt2 = interpolate.splrep(WL, RIAgIm, s=0)
RIAgImInt2 = interpolate.splev(WLx, fRIAgImInt2, der=0)

RIAu = zeros(NumWLx, dtype=complex)
epAu = zeros(NumWLx, dtype=complex)
RIAg = zeros(NumWLx, dtype=complex)
epAg = zeros(NumWLx, dtype=complex)

RIAu = RIAuReInt2 + 1j * RIAuImInt2
RIAg = RIAgReInt2 + 1j * RIAgImInt2
epAu = RIAu**2
epAg = RIAg**2
