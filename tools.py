def diffusivity_hco_in_air(T,p,nC,nH,nO):
    """
    d=diffusivity_hco_in_air(T,p,nC,nH,nO)
    diffusivity of non-ring HCO compound in air
    Use Fuiller, Schettier, and Giddings correlation
    as in Eq. 11-4.1 in Reid, Prausnitz and Sherwood

    input: T in K; p in bar;
    nC, nH and nO are the number of C, H and O atoms in the HC
    output diffusivity d in m2/s
    """
    from numpy import sqrt

    wC=12.011;wH=1.008;wO=16;wair=28.97;
    vC=16.5;vH=1.98;vO=5.48;

    vHCO=nC*vC+nH*vH+nO*vO;
    vair=20.1;
    wHCO=nC*wC+nH*wH+nO*wO;

    pinv=1./p;

    d=1e-3*T**1.75*sqrt((wHCO+wair)/(wHCO*wair))*pinv/ \
            (vHCO**.3333+vair**.3333)**2;
    return d*1e-4; #result in m2/s
    
