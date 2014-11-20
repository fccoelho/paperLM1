
# PySCeS test input file
# Stochastic Simulation Algorithm input format
# Decay Dimerazing model

Modelname: Dengue_4_serotypes
Description: Dengue model with four serotypes and cross-immunity

# Transitions
# S > I1
# S > I2
# S > I3
# S > I4
# I1 > R1
# I2 > R2
# I3 > R3
# I4 > R4
# R1 > I12
# R1 > I13
# R1 > I14
# R2 > I21
# R2 > I23
# R2 > I24
# R3 > I31
# R3 > I32
# R3 > I34
# R4 > I41
# R4 > I42
# R4 > I43
# I12 > R
# I13 > R
# I14 > R
# I21 > R
# I23 > R
# I24 > R
# I31 > R
# I32 > R
# I34 > R
# I41 > R
# I42 > R
# I43 > R



Birth:
    $pool > S
    (S+I1+I2+I3+I4+R1+R2+R3+R4+I12+I13+I14+I21+I23+I24+I31+I32+I34+I41+I42+I43+R)*mu

Death_S:
    S > $pool
    mu*S

Death_I1:
    I1 > $pool
    mu*I1

Death_I2:
    I2 > $pool
    mu*I2

Death_I3:
    I3 > $pool
    mu*I3

Death_I4:
    I4 > $pool
    mu*I4

Death_R1:
    R1 > $pool
    mu*R1

Death_R2:
    R2 > $pool
    mu*R2

Death_R3:
    R3 > $pool
    mu*R3

Death_R4:
    R4 > $pool
    mu*R4

Inf_I1:
    S > I1
    beta*S*(I1+(phi*I21)+(phi*I31)+(phi*I41))

Inf_I2:
    S > I2
    beta*S*(I2+(phi*I12)+(phi*I32)+(phi*I42))

Inf_I3:
    S > I3
    beta*S*(I4+(phi*I13)+(phi*I23)+(phi*I43))

Inf_I4:
    S > I4
    beta*S*(I4+(phi*I14)+(phi*I24)+(phi*I34))

Rec_I1:
    I1 > R1
    sigma*I1

Rec_I2:
    I2 > R2
    sigma*I2

Rec_I3:
    I3 > R3
    sigma*I3

Rec_I4:
    I4 > R4
    sigma*I4

#Secondary Infections

# R1
Inf_I12:
    R1 > I12
    R1*gamma*delta*beta*(I2+(phi*I12)+(phi*I32)+(phi*I42))

Inf_I13:
    R1 > I13
    R1*gamma*delta*beta*(I3+(phi*I13)+(phi*I23)+(phi*I43))

Inf_I14:
    R1 > I14
    R1*gamma*delta*beta*(I4+(phi*I14)+(phi*I24)+(phi*I34))

#R2
Inf_I21:
    R2 > I21
    R2*gamma*delta*beta*(I1+(phi*I21)+(phi*I31)+(phi*I41))

Inf_I23:
    R2 > I23
    R2*gamma*delta*beta*(I3+(phi*I13)+(phi*I23)+(phi*I43))

Inf_I24:
    R2 > I14
    R2*gamma*delta*beta*(I4+(phi*I14)+(phi*I24)+(phi*I34))

#R3
Inf_I31:
    R3 > I31
    R3*gamma*delta*beta*(I1+(phi*I21)+(phi*I31)+(phi*I41))

Inf_I32:
    R3 > I32
    R3*gamma*delta*beta*(I2+(phi*I12)+(phi*I32)+(phi*I42))

Inf_I34:
    R3 > I34
    R3*gamma*delta*beta*(I4+(phi*I14)+(phi*I24)+(phi*I34))

#R4
Inf_I41:
    R4 > I41
    R4*gamma*delta*beta*(I1+(phi*I21)+(phi*I31)+(phi*I41))

Inf_I42:
    R4 > I42
    R4*gamma*delta*beta*(I2+(phi*I12)+(phi*I32)+(phi*I42))

Inf_I43:
    R4 > I43
    R4*gamma*delta*beta*(I3+(phi*I13)+(phi*I23)+(phi*I43))

# Secondary Recoveries

Rec_I12:
    I12 > R
    sigma*I12

Rec_I13:
    I13 > R
    sigma*I13

Rec_I14:
    I14 > R
    sigma*I14

Rec_I21:
    I21 > R
    sigma*I21

Rec_I23:
    I23 > R
    sigma*I23

Rec_I24:
    I24 > R
    sigma*I24

Rec_I31:
    I31 > R
    sigma*I31

Rec_I32:
    I32 > R
    sigma*I32

Rec_I34:
    I34 > R
    sigma*I34

Rec_I41:
    I41 > R
    sigma*I41

Rec_I42:
    I42 > R
    sigma*I42

Rec_I43:
    I43 > R
    sigma*I43

    
#InitPar
mu = 1/70.0
beta = 400
phi = 1
sigma = 0.1
gamma = 1
delta = 1

#InitVar
S = 50000
I1 = 50
I2 = 0
I3 = 50
I4 = 50
R1 = 0
R2 = 0
R3 = 0
R4 = 0
I12 = 0
I13 = 0
I14 = 0
I21 = 0
I23 = 0
I24 = 0
I31 = 0
I32 = 0
I34 = 0
I41 = 0
I42 = 0
I43 = 0
R = 0

#Event: denv2, _TIME_ > 52 , 0 {
#I2 = 50
#}
