import glob
import copy
import pdb

def find_use_dependencies(directory="."):
   uses_dict = {}
   modules_dict = {}
   for fname in sorted(glob.glob(directory + "/*.f90")):
      f = open(fname, "r")
      flines = f.readlines()
      f.close()

      # Expanding all #include statements
      encountered_include = True
      while encountered_include:
         encountered_include = False
         for i, line in enumerate(flines):
            splitline = line.split()
            if len(splitline) > 0 and splitline[0] in ["#include", "include"]:
               f2 = open(directory + "/" + splitline[1].strip("'").strip('"'), "r")
               flines2 = f2.readlines()
               f2.close()
               encountered_include = True
               flines.pop(i)
               break
         # If we hit an include, replace it with the appropriate lines
         if encountered_include:
            flines[i:i] = flines2

      module = ""
      using = set()
      modules = set()

      for line in flines:
         splitline = line.lower().split()
         if len(splitline) < 2:
            continue
         # Keep track of modules being created
         if splitline[0] == "module" and splitline[1] != "procedure":
            modules.add(splitline[1].strip(","))
         if splitline[0] == "use":
            using.add(splitline[1].strip(","))

      routine_name = fname.split("/")[-1].replace(".f90", "")
      # if routine_name == "onetep_split":
      #    for line in flines:
      #       print line[:-1]
      #    print using
      #    print modules

      uses_dict[routine_name] = using
      modules_dict[routine_name] = modules
   return uses_dict, modules_dict

def print_dependencies(routine, modules, max_width=100):
   next_line = routine + ".out:"
   indent = 3
   for string in [routine + ".f90"] + [m + ".o" for m in modules]:
      if len(next_line + string) > max_width - 3:
         print next_line + " \\"
         next_line = " "*(indent)
      next_line += " " + string
   print next_line

def merge_two_dicts(x, y):
    z = x.copy()   # start with x's keys and values
    z.update(y)    # modifies z with y's keys and values & returns None
    return z

def obtain_all_dependencies(module, dependencies, direct_dependencies, depth=0):
   # module - the module whose dependencies we're interested in
   # dependencies - a set of dependencies (both direct and indirect) of module
   # direct_dependencies - the dictionary of direct_dependencies
   # extra_dependencies - returns the extra dependencies not already in dependencies
   # print "  "*depth + module
   # print sorted(dependencies)
   # if depth > 6:
   #    quit()
   if module in direct_dependencies.keys():
      for dep_mod in direct_dependencies[module]:
         if dep_mod in dependencies:
            continue
         dependencies.update(obtain_all_dependencies(dep_mod, dependencies, direct_dependencies, depth+1))
   dependencies.add(module)
   return dependencies
   
if __name__ == "__main__":

   uses_utils, _ = find_use_dependencies("../lib/utils")
   uses_ed_solver, _ = find_use_dependencies("../ed_solver")
   uses_scripts, scripts_modules = find_use_dependencies()
   direct_dependencies = merge_two_dicts(uses_ed_solver, uses_scripts)
   direct_dependencies = merge_two_dicts(direct_dependencies, uses_utils)

   for script, modules in sorted(uses_scripts.iteritems()):
      all_uses = obtain_all_dependencies(script, set(), direct_dependencies)

      # Do not include use statements for itself, mpi, or any modules defined within the file itself
      for module in list(scripts_modules[script]) + [script, "mpi"]:
         try:
            all_uses.remove(module)
         except:
            pass
      print_dependencies(script, sorted(list(all_uses)))
      print ""
