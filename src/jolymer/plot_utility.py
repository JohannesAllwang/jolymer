# -*- coding: utf-8 -*-
"""
Created on Tue Dec 29 10:21:01 2020

@author: xcill
"""

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mc
import colorsys

    
def colormap(name, start, stop, num):
    cm = plt.get_cmap(name)
    listed = cm(np.linspace(start, stop, num))
    out = iter(listed)
    return out
def cm_for_l(name, l, start=0.1, stop=0.9):
    num = len(l)
    return colormap(name, start, stop, num)

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
