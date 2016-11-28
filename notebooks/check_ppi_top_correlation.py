def main():
  load_ppi()

def load_ppi():
  filename = '../PPI/PPI.sig'

  f = open(filename, 'r')
  lines = f.readlines()
  f.close()

  print(len(lines))

  for i in range(10):

    if i > 0:
      inst_line = lines[i]

      inst_line = inst_line.strip().split()

      source = inst_line[0]

      target = inst_line[5]


      print('\n')
      print(inst_line)
      print(source)
      print(target)


main()