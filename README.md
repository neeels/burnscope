burnscope
=========

Mesmerizing colorful animation from an underdamped blur algorithm.
(c) 2014 Neels Hofmeyr <neels@hofmeyr.de>

Published under the GNU General Public License v3.
Burnscope produces a mesmerizing animation that I discovered by accident when I
was a teenager. I've recreated it in memories of old times. It repeatedly
applies a simple underdamped blur algorithm to a seed image, allowing the color
values to wrap when overflowing. If you can explain how this staggering
everchanging complexity can spring from such a simple algorithm and just one
pixel as seed, please send me an email ;)

Usage example:
  burnscope -g 320x200 -m 2 -p 70
  burnscope3 -g 320x200 -m 2 -p 70

burnscope uses an indexed palette and runs a single burnscope to cycle through
the color palette.

burnscope3 uses three independent burnscopes to feed the color channels r, g,
and b.

To run, you need to have the SDL development packages installed. Simply run
'make' to build, and run ./burnscope or ./burnscope3 to launch an animation.
The -h option prints a commandline help.

Many enhancements come to mind. If you enjoy playing with this, I would be glad
to hear about it! Patches welcome!
