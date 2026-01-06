
# gui/regals/plots/uv.py
def plot_uv_alignment(cm_saxs, ax):
    ax.clear()
    cm_saxs.plot_alignment(ax=ax)
    ax.set_title("UV–SAXS alignment")
    ax.figure.canvas.draw_idle()


def plot_interpolated_uv(cm_saxs, ax):
    ax.clear()
    uv = cm_saxs.interpolate_uv()
    ax.plot(uv["x"], uv["y"])
    ax.set_title("Interpolated UV")
    ax.figure.canvas.draw_idle()
