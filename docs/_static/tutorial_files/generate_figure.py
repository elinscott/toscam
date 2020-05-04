#!/usr/bin/pymol

# Script for automatically generating figure of FeCO3, with nice colours
from pymol import cmd
from pymol.vfont import plain

def select_by_id_list(sele_name, id_list):
   sele_string = "index " + " | index ".join(map(str,id_list))
   cmd.select(sele_name, sele_string)

def add_label(label, position, label_name = None):
   if label_name is None:
      label_name = label + "_label"
   text = [COLOR, 1.0, 1.0, 0.0,]
   w = 0.06
   axes = 4*w*rotation
   position -= w*camera_rightwards
   position -= w*camera_upwards
   cyl_text(text,plain,position,label,0.5*w,axes=axes)
   cmd.load_cgo(text, label_name)

# Loading
cmd.load("FeCO3.xyz")

cmd.set_view ([\
     0.454534203,   -0.583148003,    0.673300803,\
     0.870026767,    0.128641263,   -0.475922495,\
     0.190919206,    0.802114487,    0.565826774,\
    -0.000015825,   -0.000005098,   -8.543414116,\
     5.045799255,    5.481069088,    5.299679756,\
  -210.128402710,  227.214721680,  -20.000000000])

# Figure settings
cmd.hide("lines")
cmd.set("transparency", 0.2)
cmd.set("spec_reflect", 0)
cmd.set("spec_power", 1500)

# Coloring
cmd.bg_color("white")
seles = ["carbons", "oxygens"]
eles = ["C", "O"]
# see https://pymolwiki.org/index.php/Color_Values#Simple_named_colours
colors = ["lutetium", "oxygen"]
for sele, ele, color in zip(seles, eles, colors):
   cmd.select(sele, "elem " + ele)
   cmd.color(color, sele)

# Display mode
cmd.show("sticks")
cmd.set("stick_radius", 0.14)
cmd.set("sphere_scale", 0.25)
cmd.show("spheres")

cmd.ray(600, 500)
cmd.png("feco3.png")
