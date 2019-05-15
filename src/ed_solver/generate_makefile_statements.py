import glob

if __name__ == "__main__":
   out = []
   for fname in sorted(glob.glob("*.f90")):
      f = open(fname, "r")
      flines = f.readlines()
      f.close()

      module = ""
      using = set()

      for line in flines:
         splitline = line.lower().split()
         if len(splitline) < 2:
            continue
         if splitline[0] == "module" and splitline[1] != "procedure":
            module = splitline[1]
         if splitline[0] == "use":
            using.add(splitline[1].strip(","))

      moduleList = ""
      if len(using) > 0:
         moduleList = ".o ".join(sorted(list(using))) + ".o"

      out.append("{0}.o: {0}.f90 ".format(module) + moduleList)
   

   # Print routines:
   for o in sorted(out):
      print o.split(":")[0].replace(".mod", ".f90")

   # Print dependencies:
   for o in sorted(out):
      print o
      print ""
