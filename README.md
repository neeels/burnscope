burnscope
=========

Mesmerizing colorful animation from an underdamped blur algorithm.

(c) 2014,2016 Neels Hofmeyr <neels@hofmeyr.de>

Published under the GNU General Public License v3.

![screenshot1](http://kleinekatze.de/no4uF3Pa/burnscope_screenshot.png)
![screenshot2](http://kleinekatze.de/no4uF3Pa/burnscope_screenshot2.png)
![screenshot3](http://kleinekatze.de/no4uF3Pa/burnscope_screenshot3.png)

Burnscope produces a mesmerizing animation that I discovered by accident when I
was a teenager. I've recreated it in memories of old times. It repeatedly
applies a simple underdamped blur algorithm to a seed image, allowing the color
values to wrap when overflowing. If you can explain how this staggering
everchanging complexity can spring from such a simple algorithm and just one
pixel as seed, please send me an email ;)

burnscope\_fft is the latest greatest: it's fast and has numerous features the
other two (burnscope and burnscope3) don't have.

With burnscope\_fft, you can *attach a game controller* and influence the
animation. You can save to raw video file with music playing along (to make a
music video). See the -h option.

    ./burnscope_fft -h

Usage examples:

    sudo apt-get install gcc libsdl1.2-dev libpng12-dev libfftw3-dev libsndfile1-dev
    cd burnscope
    make
    ./burnscope_fft -g 320x200 -m 2 -f 25
    ./burnscope3 -g 320x200 -m 2 -p 70

To find out all features, you'll have to read the source code:

* keyboard shortcuts
* PNG images dropped into burnscope\_fft
* ...

burnscope3 uses three independent burnscopes to feed the color channels r, g,
and b. It's not as great as I had expected ;)

Many enhancements come to mind. If you enjoy playing with this, I would be glad
to hear about it!

Patches welcome!
