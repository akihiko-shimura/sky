# -*- coding: utf-8 -*-
"""
Created on Wed Jun 01 02:56:37 2016

@author: A.Shimura

fs_kHzで取得した実験データをPythonで便利に扱うためのデータ構造のclass

Wdat  :  スペクトルデータ(in Wavelength[nm])。
Fdat  :  スペクトルデータ(in Frequency[THz])。Wdatから変換される。
Rdat  :  ラマンシフトデータ(in Raman shift[cm-1])。ラマン励起波長を与えてWdatから変換される。
Tdat  :  時間データ（in Time[ps]）。Btach WdataをTconversionしたデータ。

"""

#import linecache
import numpy
import os
from scipy.interpolate import interp1d

from . import helper

__all__ = ["Wdat", "Fdat", "Rdat", "Tdat"]

# データファイルの形式の定義
SIG_TYPE_LINE_NUMBER = 8  # sky data fileでsignal data typeが書かれている行番号


class Wdat(object):
    """波長軸のスペクトルデータ"""
    __slot__ = ["delay", "x", "y", "name", "index"]

    def __init__(self):
        self.delay = None # delay time(s)
        self.x = numpy.array([]) # wavelength array
        self.y = numpy.array([]) # signal array
        self.name = 'wdat'
        self.index = 0

    def load(self, fpath):
        """既存のWdata fileをloadする"""
        """
        headerからdelay times情報を得る。
        Format: wdat.delay = [0.00, 0.50, 0.00, 0.00]
        """
        with open(fpath, 'r') as f:
            for _line in f.readlines():
                """#EXで始まる行"""
                if _line.startswith('#EX'):
                    self.delay = list(map(float, _line.strip().split()[1:]))

#        _line_of_delays = getline(fpath, DELAY_NUMBER) # 2nd line
#        try:
#            self.delay = map(float, _line_of_delays.strip().split()[1:])
#        except:
#            self.delay = [0,0,0,0]

        """file name and index
        file name "DA1050_DA_0001"の0001の部分"""
        self.name = os.path.basename(fpath)
        self.index = int(self.name.split('_')[-1]) # last element, int or str

        """x, y array"""
        _xy = numpy.loadtxt(fpath)
        self.x = _xy[:,0] # 0th column
        self.y = _xy[:,1] # 1st column

        """signalのtype（SG1, SG2, DA etc.）の情報"""
#        _line_sig_type = linecache.getline(fpath, SIG_TYPE_LINE_NUMBER) #8th line
#        self.sig_type = _line_sig_type[1:].strip()

        return self

#    @property
#    def _yFunc(self, kind='cubic'):
#        """
#        Wdat.yを補完した関数を返す。
#        注意：引数はそれが単調増加でなければならない。
#        """
#        return scipy.interpolate.interp1d(self.x, self.y, kind,
#                        bounds_error=False, fill_value=0)


_nm2THz = lambda x : 299792.45800 / x # c/x(nm)
_nm2eV  = lambda x : 1240.0 / x

class Fdat(object):
    """
    Wdat objectを対応する等間隔周波数軸のobjectに変換する。単位はnmからTHzに変換される。
    逆変換（THzからnm）も同じ。
    __USAGE
    wd = Wdat(0.5, x) # 親object
    fd = Fdat(wd)
    __FUTURE
    親クラスの属性へのアクセスをsuper(Fdat, self)としてもいい。そうすればWdatとあからさま
    に指定しなくてよく，より汎用化される。
    """
    __slot__ = ["delay", "x", "y", "name", "index"]

    def __init__(self, wdat, unit='THz'):
        self.delay = wdat.delay

        if unit == 'THz':
            x_thz_equistep = numpy.linspace(_nm2THz(wdat.x[-1]),
                                            _nm2THz(wdat.x[0]),
                                            len(wdat.x))
            self.x = x_thz_equistep
            """wdat.yを補完した関数の引数とする非等間隔x"""
            x_nm_new = _nm2THz(x_thz_equistep)

        if unit == 'eV':
            x_eV_equistep = numpy.linspace(_nm2eV(wdat.x[-1]),
                                           _nm2eV(wdat.x[0]),
                                           len(wdat.x))
            self.x = x_eV_equistep
            """wdat.yを補完した関数の引数とする非等間隔x"""
            x_nm_new = _nm2eV(x_eV_equistep)

        """逆変換した_newx(nm)がfloatの誤差でinterpolate rangeの外に出る問題に対処する"""
        # old way
#        x_nm_new[0]  = wdat.x[0]
#        x_nm_new[-1] = wdat.x[-1]
#        yfunc = interp1d(wdat.x, wdat.y, kind=kind,
#                         bounds_error=False, fill_value=0)
        # new way
        yfunc = interp1d(wdat.x, wdat.y, kind='linear',
                         bounds_error=False, fill_value='extrapolate')

        self.y = yfunc(x_nm_new)
#        self.y = wdat._yFunc()(_x_nm_new)

        self.name = wdat.name
        self.index = int(wdat.index)


class Rdat(object):
    """
    波長軸のスペクトルデータから変換されるRaman shift軸のラマンスペクトルデータ

    Parameter
    ---------
    x_Rpump  :  Raman pump波長 [nm]

    Wdat objectを対応する等間隔Ramans shift軸のobjectに変換する。
    xの単位はcm-1になる。x > 0 はStokes, x < 0 はanti-Stokes散乱。

    __USAGE
    wd = Wdat(0.2, x) # ある波長領域スペクトルを作る。
    rd = Rdat(wd, x_Rpump) # Raman shift表示した等間隔スペクトル
    """
    __slot__ = ["delay", "x", "y", "name", "index", "x_Rpump"]

    def __init__(self, wdat, x_Rpump):
        self.x_Rpump = x_Rpump

        self.delay = wdat.delay
        q_equistep = numpy.linspace(helper.raman_shift(wdat.x[-1], x_Rpump),
                                 helper.raman_shift(wdat.x[0], x_Rpump),
                                 len(wdat.x)) #cm-1
        self.x = q_equistep
        """wdat.yを補完した関数の引数とする非等間隔x"""
        x_nm_new = helper.inverse_raman_shift(q_equistep, x_Rpump) #nm
        """逆変換した_newx(nm)がfloatの誤差でinterpolate rangeの外に出る問題に対処する"""
#        x_nm_new[0]  = wdat.x[0]
#        x_nm_new[-1] = wdat.x[-1]
#        yfunc = interp1d(wdat.x, wdat.y, kind=kind,
#                         bounds_error=False, fill_value=0)
        yfunc = interp1d(wdat.x, wdat.y, kind='linear',
                         bounds_error=False, fill_value='extrapolate')
        self.y = yfunc(x_nm_new)
#        self.y = wdat._yFunc(kind)(x_nm_new)

        self.name = wdat.name
        self.index = int(wdat.index)

#==============================================================================

class Tdat(object):
    """時間軸に変換したダイナミクスデータ"""
    __slot__ = ["wl", "t", "y", "name", "index"]
    def __init__(self):
        self.wl = None # wavelength(s)
        self.t = numpy.array([]) # delay time array
        self.y = numpy.array([]) # signal array
        self.name = 'tdat'
        self.index = 0

    def load(self, fpath):
        with open(fpath, 'r') as f:
            for _line in f.readlines():
                """headerからwavelengths情報を得る"""
                if _line.startswith('#EX'):
                    self.wl = list(map(float, _line.strip().split()[1:]))

        """x, y array"""
        _xy = numpy.loadtxt(fpath)
        self.t = _xy[:,0] # 0th column
        self.y = _xy[:,1] # 1st column

        """file name and index (file name: DA1050T.001 の 001 の部分)"""
        self.name = os.path.basename(fpath)

        try:
            self.index = int(self.name.split('.')[-1]) # last element

        except ValueError:
            # for mean Tdata, e.g., mean_620-nm_027.dat --> 027
            self.index = int(self.name.split('.')[-2])

    #    """
    #    signalのtype（SG1, SG2, DA etc.）の情報
    #    """
    #    _line7 = linecache.getline(fpath, 3) # 7th line
    #    tdat.sig_type = _line7[1:4].strip()

        return self
