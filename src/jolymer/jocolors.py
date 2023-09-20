# https://colordesigner.io/#0D6986-0000E0-D2998E-86660D-5B1A07
import matplotlib.colors as mc

#0D6986
#0000E0
#D2998E
#86660D
#5B1A07

lllCyan = '#9ce0f5'
llCyan = '#8ec7d2'
Cyan = '#0d6986'
dCyan = '#07485b'
ddCyan = '#042028'

lllRed = '#f9cbbd'
llRed = '#ef7d5a'
Red = '#862a0d'
dRed = '#431507'
ddRed = '#1b0803'
Blue1 = '#0d2d86'
Blue2 = '#2a0d86'
llBlue = '#7d5aef'
Blue = '#2a0d86'
dBlue = '#150743'
# llBlue = '#998ed2'
# Blue = '#3f13c8'
# dBlue = '#1a075b'
llRust = '#d2998e'
Rust = '#ed4919'
dRust = '#5b1a07'
llPurple = '#bb8ed2'
Purple = '#660d86'
dPurple = '#44075b'
Magenta = '#860d69'
llGreen = '#a5d28e'
Green = '#2d860d'
dGreen = '#1e5b07'
GreenYellow = '#69860d'
llBrown = '#d2bb8e'
Brown = '#86660d' # split complementary
dBrown = '#5b4407'


Yellow = '#e0e000'
Orange = '#c57e00'
Wine = '#860d2d'
# Red = '#d50000'

White = '#FFFFFF'
Black = '#000000'
Gray1 = '#d0d0d0'
Gray2 = '#f0f0f0'

Complementary = Red
Analogous = [Blue1, Blue2]
SplitComplementary = [Wine, Brown]
Triad = [Magenta, GreenYellow]
Square = [Purple, Red, Green]
Rectangle = [Blue2, Red, GreenYellow]

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    import numpy as np
    primary = Cyan
    tri1 = f'#{primary[5:7]}{primary[1:3]}{primary[3:5]}'
    tri2 = f'#{primary[3:5]}{primary[5:7]}{primary[1:3]}'
    tri2 = '#860d2d'
    tri3 = f'#{primary[1:3]}{primary[5:7]}{primary[3:5]}'
    tri3 = f'#86660d'
    testmap = mc.LinearSegmentedColormap.from_list('test', [primary, tri2, tri3, primary])
    testmap = mc.LinearSegmentedColormap.from_list('test', [Red, 'white', Blue])
    colors = testmap(np.linspace(0, 1, num=25))
    # colors = [Cyan, '#862a0d', Purple, Green]
    x = np.linspace(0, 10, num=15)
    y = np.linspace(0, 10, num=15)
    for color in colors:
        plt.plot(x, y, '-o', color=color)
        y += 1
    plt.show()
