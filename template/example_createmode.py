from ramanbond.adfnormalmodes import *


print('Use the normalmode analysis feature of the ramanbond package to create the input files of ADF aoresponse.')
print('A vibration calculation output file containing information of normal modes and a template file defining the structure of new files are needed.')
print('The example systems is CO. The vibration calculation output is freq_CO.out and the template file is aoresponse_template.run. ')
vib_CO = adfnormalmodes('freq_CO.out')

print('Define the vibrational frequency range:')
freqmin = 1000
freqmax = 3000

print('The studied frequencies range from {0: .3f} to {1: .3f}'.format(freqmin, freqmax))

for freq in vib_CO.freq:
  if freq > freqmin and freq < freqmax:
    vib_CO.create_modefile(mode = freq, num_step = 1, template = 'aoresponse_template.run')

print('The input files of ADF aoresponse are created.')

