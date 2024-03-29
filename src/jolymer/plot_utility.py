# -*- coding: utf-8 -*-
"""
Created on Tue Dec 29 10:21:01 2020

@author: xcill
"""

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mc
import colorsys


# For plotting:
tum_purple = '#69085a'
tum_dblue = '#0f1b5f'
tum_cyan = '#00778a'
tum_dgreen = '#007c30'
tum_green = '#679a1d'
tum_yellow = '#ffdc00'
tum_dyellow = '#f9ba00'
tum_orange = '#d64c13'
tum_red = '#c4071b'
tum_dred = '#9c0d16'
# primary colors
tum_p1 = '#0065bd'
tum_blue = '#0065bd'
# secondary colors:
tum_dblue = '#005293'
tum_ddblue = '#003359'
tum_s3 = '#333333'
tum_s4 = '#808080'
tum_s5 = '#cccccc'
# accent colors
tum_a1 = '#dad7cb'
tum_a2 = '#e37222'
tum_aorange = '#e37222'
tum_a3 = '#a2ad00'
tum_a4 = '#98c6ea'
tum_a5 = '#64a0c8'

tum_colors = [
    tum_dblue,
    tum_cyan,
    tum_dgreen,
    tum_green,
    # tum_yellow,
    tum_dyellow,
    tum_orange,
    tum_red,
    tum_dred,
    tum_purple]


def loglog(ax):
    """necessary?"""
    ax.set_yscale('log')
    ax.set_xscale('log')


def alterhex(color):
    """
    Make a color class!
    """
    return color.replace('#', '0x')


def colormap(name, start, stop, num):
    """
    better use the cycler?
    """
    cm = plt.get_cmap(name)
    listed = cm(np.linspace(start, stop, num))
    out = iter(listed)
    return out


def cm_for_l(name, list, start=0.1, stop=0.9):
    """
    remove this?
    """
    num = len(list)
    return colormap(name, start, stop, num)


# top = tum_blues
# bottom = cm.get_cmap('Blues', 128)
# newcolors = np.vstack((top(np.linspace(0, 1, 128)),
                       # bottom(np.linspace(0, 1, 128))))

tum_dblue_white_dred = mc.LinearSegmentedColormap.from_list("tum_blue_to_red", [tum_dred, 'white', tum_dblue])
tum_blue_white_red = mc.LinearSegmentedColormap.from_list("tum_blue_to_red", [tum_red, 'white', tum_blue])


def colorgradient(name, colors, num):
    lscm = mc.LinearSegmentedColormap.from_list(name, colors)
    listed = lscm(np.linspace(0, 1, num))
    out = iter(listed)
    return out




def n_subplots(n, subplot_size=(5, 4)):
    sw, sh = subplot_size
    ncols = min(n, 2)
    nrows = int(np.ceil(n/2))
    fig, axes = plt.subplots(ncols = ncols, nrows=nrows,
                           squeeze= False, figsize=(sw*ncols, sh*nrows))
    return fig, axes
def lighten_color(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])


tum_lblue = lighten_color(tum_blue, amount=0.7)
tum_lpurple = lighten_color(tum_purple, amount=0.7)
tum_lcyan = lighten_color(tum_cyan, amount=0.7)
tum_lgreen = lighten_color(tum_green, amount=0.7)
tum_lyellow = lighten_color(tum_yellow, amount=0.7)
tum_ldyellow = lighten_color(tum_dyellow, amount=0.7)
tum_lorange = lighten_color(tum_orange, amount=0.7)
tum_lred = lighten_color(tum_red, amount=0.7)
tum_ldred = lighten_color(tum_dred, amount=0.7)

tum_llblue = lighten_color(tum_blue, amount=0.3)
tum_llpurple = lighten_color(tum_purple, amount=0.3)
tum_llcyan = lighten_color(tum_cyan, amount=0.3)
tum_llgreen = lighten_color(tum_green, amount=0.3)
tum_llyellow = lighten_color(tum_yellow, amount=0.3)
tum_lldyellow = lighten_color(tum_dyellow, amount=0.3)
tum_llorange = lighten_color(tum_orange, amount=0.3)
tum_llred = lighten_color(tum_red, amount=0.3)
tum_lldred = lighten_color(tum_dred, amount=0.3)

tum_blues = mc.ListedColormap(
        [tum_dblue, lighten_color(tum_dblue, 0.8), lighten_color(tum_dblue, 0.5), lighten_color(tum_dblue, 0.2)]
        )
tum_blue_to_red = mc.LinearSegmentedColormap.from_list("tum_blue_to_red", [tum_dred, lighten_color(tum_dred, 0.3), lighten_color(tum_dblue, 0.3), tum_dblue])


class data_linewidth_plot():

    """
    this class is copy pasted from stack overflow

    """
    def __init__(self, x, y, **kwargs):
        self.ax = kwargs.pop("ax", plt.gca())
        self.fig = self.ax.get_figure()
        self.lw_data = kwargs.pop("linewidth", 1)
        self.lw = 1
        self.fig.canvas.draw()

        self.ppd = 72./self.fig.dpi
        self.trans = self.ax.transData.transform
        self.linehandle, = self.ax.plot([],[],**kwargs)
        if "label" in kwargs: kwargs.pop("label")
        self.line, = self.ax.plot(x, y, **kwargs)
        self.line.set_color(self.linehandle.get_color())
        self._resize()
        self.cid = self.fig.canvas.mpl_connect('draw_event', self._resize)

    def _resize(self, event=None):
        lw =  ((self.trans((1, self.lw_data))-self.trans((0, 0)))*self.ppd)[1]
        if lw != self.lw:
            self.line.set_linewidth(lw)
            self.lw = lw
            self._redraw_later()

    def _redraw_later(self):
        self.timer = self.fig.canvas.new_timer(interval=10)
        self.timer.single_shot = True
        self.timer.add_callback(lambda : self.fig.canvas.draw_idle())
        self.timer.start()
