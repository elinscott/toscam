import ipdb

if __name__ == '__main__':
   keywords = {}

   files_with_keywords = ["scripts/onetep_dmft.f90",
                          "scripts/dmft_one_iteration_split_module.h",
                          "scripts/onetep_dmft_one_iter.f90",
                          "scripts/onetep_dmft_sc.f90",
                          "ed_solver/solver_init.h",
                          "ed_solver/dmft_ed_solver_variables.h"]

   for fname in files_with_keywords:
      with open('../src/' + fname, 'r') as f:
         flines = f.readlines()

      for i, line in enumerate(flines):
         line = line.strip()

         # Finding calls to putel_in_namelist
         if 'putel_in_namelist' not in line.lower() or line[:4].lower() != 'call':
            continue

         # Dealing with putel_in_namelist calls running over multiple lines
         while line[-1] == '&':
            i += 1
            line += flines[i]
            line = line.strip()

         line = line.replace('&','')
         line = line.rstrip(')')

         keyword, default, definition = [s.strip(' "' + "'") for s in line.split(',')[2:5]]

         # Don't include keywords whose defaults are not hard-coded
         if 'trim(' in default.lower():
            continue

         # Tidying keywords
         keyword = keyword.lower()

         # Tidying defaults
         key_type = 'string'
         if default.lower() in ['.true.', '.false.']:
            default = bool(default.strip('.'))
            key_type = 'boolean'
         elif any([c.isdigit() for c in default]) and not any([c.isalpha() and c not in ['d', 'e'] for c in default]):
            if all([c.isdigit() or c == '-' for c in default]):
               default = int(default)
               key_type = 'integer'
            else:
               default = float(default.replace('d','e'))
               key_type = 'float'

         if key_type == 'string':
            default = f"'{default}'"

         # Tidying definitions
         while '  ' in definition:
            definition = definition.replace('  ', ' ')
         if definition[:10].lower() == 'definition':
            definition = definition[10:]
            definition = definition.lstrip(': ')

         if keyword in keywords:
            print(f'WARNING: multiple definitions found for {keyword}')

         keywords[keyword] = [keyword, definition, default, key_type]

   header = ['<table id="keywordTable" style="width:100%; text-align:left">',
             '   <tr>',
             '      <th>Keyword</th>',
             '      <th>Description</th>',
             '      <th>Default</th>',
             '      <th>Type</th>',
             '   </tr>']

   with open('_static/keywords.html', 'w') as f:
      f.write('\n'.join(header))
      for key in sorted(keywords.keys()):
         keyword, definition, default, key_type = keywords[key]
         f.write(f'\n   <tr>')
         f.write(f'\n      <td><code>{keyword}</code></td>')
         f.write(f'\n      <td>{definition}</td>')
         f.write(f'\n      <td><code>{default}</code></td>')
         f.write(f'\n      <td>{key_type}</td>')
         f.write(f'\n   </tr>')
      f.write('\n</table>')
