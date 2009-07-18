from reactor import *

a = Reactor('chem.cti')
print 'the mechanism has %d species and %d reactions' % \
(len(ctml._speciesnames),len(ctml._reactions))
#print ctml._speciesnames

a.setVolFrac({'O2(1)':0.01,
    'n-undecane(2)':0.04,
    'n-tridecane(3)':0.19,
    'SPC(4)':0.11,
    'n-hexadecane(5)':0.25,
    'SPC(6)':0.12,
    'n-nonadecane(7)':0.18,
    'n-heneicosane(8)':0.10,
    'HO2J(9)':0.0,
    'C16H25J(42)':0.0,
    'C16H25O2J(117)':0.0,
    'C16H26O2(119)':0.0,
    'H2O2(101)':0.0,
    'C16H24(104)':0.0,
    '(DEPOSIT)':0.0
    })

print a.speciesIndex
print a.volFrac
print 'the volFrac of O2(1) is %f' % a.volFrac[a.speciesIndex['O2(1)']]


# g/mol
a.setMolWeight({'O2(1)':32,
    'n-undecane(2)':156.31,
    'n-tridecane(3)':184.36,
    'SPC(4)':142.20,
    'n-hexadecane(5)':226.44,
    'SPC(6)':218.378,
    'n-nonadecane(7)':268.521,
    'n-heneicosane(8)':296.574,
    'HO2J(9)':33,
    'C16H25J(42)':217.37,
    'C16H25O2J(117)':249.37,
    'C16H26O2(119)':250.37,
    'H2O2(101)':34,
    'C16H24(104)':216.36,
    '(DEPOSIT)':260.54
    })

print 'the MW of O2(1) is %f' % a.molWeight[a.speciesIndex['O2(1)']]
# g/cm^3, value assumed
a.setMassDensity({'O2(1)':0.00143,
    'n-undecane(2)':0.74,
    'n-tridecane(3)':0.7635,
    'SPC(4)':0.74,
    'n-hexadecane(5)':0.25,
    'SPC(6)':0.12,
    'n-nonadecane(7)':0.18,
    'n-heneicosane(8)':0.10,
    'HO2J(9)':0.76,
    'C16H25J(42)':0.76,
    'C16H25O2J(117)':0.87,
    'C16H26O2(119)':0.98,
    'H2O2(101)':0.87,
    'C16H24(104)':0.89,
    '(DEPOSIT)':0.90
    })
print 'the Mass Density of O2(1) is %f' % \
        a.massDens[a.speciesIndex['O2(1)']]

a.getConc()
print 'the concentration of n-undecane(2) is %f' % \
        a.concs[a.speciesIndex['n-undecane(2)']]

# need to test the ode
a.advance(1.01)
print a.concs
print 'the concentration of n-undecane(2) is %f' % \
        a.concs[0,a.speciesIndex['n-undecane(2)']]
