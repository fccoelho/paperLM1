
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
    (S+I1+I2+I3+I4+R1+R2+R3+R4+I12+I13+I14+I21+I23+I24+I31+I32+I34+I41+I42+I43+R12+R13+R14+R23+R24+R34+I231+I241+I341+I132+I142+I342+I123+I143+I243+I124+I134+I234+R123+R124+R134+R234+I1234+I1243+I1342+I2341+R1234)*mu

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
    
Death_I12:
    I12 > $pool
    mu*I12

Death_I13:
    I13 > $pool
    mu*I13

Death_I14:
    I14 > $pool
    mu*I14
    
Death_I21:
    I21 > $pool
    mu*I21
    
Death_I23:
    I23 > $pool
    mu*I23

Death_I24:
    I24 > $pool
    mu*I24

Death_I31:
    I31 > $pool
    mu*I31

Death_I32:
    I32 > $pool
    mu*I32

Death_I34:
    I34 > $pool
    mu*I34

Death_I41:
    I41 > $pool
    mu*I41

Death_I42:
    I42 > $pool
    mu*I42

Death_I43:
    I43 > $pool
    mu*I43

Death_R12:
    R12 > $pool
    mu*R12
    
Death_R13:
    R13 > $pool
    mu*R13
    
Death_R14:
    R14 > $pool
    mu*R14
    
Death_R23:
    R23 > $pool
    mu*R23
    
Death_R24:
    R24 > $pool
    mu*R24
    
Death_R34:
    R34 > $pool
    mu*R34
    
Death_I231:
    I231 > $pool
    mu*I231

Death_I241:
    I241 > $pool
    mu*I241

Death_I341:
    I341 > $pool
    mu*I341
    
Death_I132:
    I132 > $pool
    mu*I132
    
Death_I142:
    I142 > $pool
    mu*I142

Death_I342:
    I342 > $pool
    mu*I342

Death_I123:
    I123 > $pool
    mu*I123

Death_I143:
    I143 > $pool
    mu*I143

Death_I243:
    I243 > $pool
    mu*I243

Death_I124:
    I124 > $pool
    mu*I124

Death_I134:
    I134 > $pool
    mu*I134

Death_I234:
    I234 > $pool
    mu*I234

Death_R123:
    R123 > $pool
    mu*R123

Death_R124:
    R124 > $pool
    mu*R124

Death_R134:
    R134 > $pool
    mu*R134

Death_R234:
    R234 > $pool
    mu*R234

Death_I1234:
    I1234 > $pool
    mu*I1234

Death_I1243:
    I1243 > $pool
    mu*I1243

Death_I1342:
    I1342 > $pool
    mu*I1342

Death_I2341:
    I2341 > $pool
    mu*I2341

Death_R1234:
    R1234 > $pool
    mu*R1234
    
# Primary Infections
    
Inf_I1:
    S > I1
    beta*S*(I1+(phi*I21)+(phi*I31)+(phi*I41)+(phi*I231)+(phi*I241)+(phi*I341)+(phi*I2341))

Inf_I2:
    S > I2
    beta*S*(I2+(phi*I12)+(phi*I32)+(phi*I42)+(phi*I132)+(phi*I142)+(phi*I342)+(phi*I1342))

Inf_I3:
    S > I3
    beta*S*(I3+(phi*I13)+(phi*I23)+(phi*I43)+(phi*I123)+(phi*I143)+(phi*I243)+(phi*I1243))

Inf_I4:
    S > I4
    beta*S*(I4+(phi*I14)+(phi*I24)+(phi*I34)+(phi*I124)+(phi*I134)+(phi*I234)+(phi*I1234))
    
# Primary Recoveries

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
    R1*gamma*delta*beta*(I2+(phi*I12)+(phi*I32)+(phi*I42)+(phi*I132)+(phi*I142)+(phi*I342)+(phi*I1342))

Inf_I13:
    R1 > I13
    R1*gamma*delta*beta*(I3+(phi*I13)+(phi*I23)+(phi*I43)+(phi*I123)+(phi*I143)+(phi*I243)+(phi*I1243))

Inf_I14:
    R1 > I14
    R1*gamma*delta*beta*(I4+(phi*I14)+(phi*I24)+(phi*I34)+(phi*I124)+(phi*I134)+(phi*I234)+(phi*I1234))

#R2
Inf_I21:
    R2 > I21
    R2*gamma*delta*beta*(I1+(phi*I21)+(phi*I31)+(phi*I41)+(phi*I231)+(phi*I241)+(phi*I341)+(phi*I2341))

Inf_I23:
    R2 > I23
    R2*gamma*delta*beta*(I3+(phi*I13)+(phi*I23)+(phi*I43)+(phi*I123)+(phi*I143)+(phi*I243)+(phi*I1243))

Inf_I24:
    R2 > I24
    R2*gamma*delta*beta*(I4+(phi*I14)+(phi*I24)+(phi*I34)+(phi*I124)+(phi*I134)+(phi*I234)+(phi*I1234))

#R3
Inf_I31:
    R3 > I31
    R3*gamma*delta*beta*(I1+(phi*I21)+(phi*I31)+(phi*I41)+(phi*I231)+(phi*I241)+(phi*I341)+(phi*I2341))

Inf_I32:
    R3 > I32
    R3*gamma*delta*beta*(I2+(phi*I12)+(phi*I32)+(phi*I42)+(phi*I132)+(phi*I142)+(phi*I342)+(phi*I1342))

Inf_I34:
    R3 > I34
    R3*gamma*delta*beta*(I4+(phi*I14)+(phi*I24)+(phi*I34)+(phi*I124)+(phi*I134)+(phi*I234)+(phi*I1234))

#R4
Inf_I41:
    R4 > I41
    R4*gamma*delta*beta*(I1+(phi*I21)+(phi*I31)+(phi*I41)+(phi*I231)+(phi*I241)+(phi*I341)+(phi*I2341))

Inf_I42:
    R4 > I42
    R4*gamma*delta*beta*(I2+(phi*I12)+(phi*I32)+(phi*I42)+(phi*I132)+(phi*I142)+(phi*I342)+(phi*I1342))

Inf_I43:
    R4 > I43
    R4*gamma*delta*beta*(I3+(phi*I13)+(phi*I23)+(phi*I43)+(phi*I123)+(phi*I143)+(phi*I243)+(phi*I1243))

# Secondary Recoveries

Rec_I12:
    I12 > R12
    sigma*I12

Rec_I13:
    I13 > R13
    sigma*I13

Rec_I14:
    I14 > R14
    sigma*I14

Rec_I21:
    I21 > R12
    sigma*I21

Rec_I23:
    I23 > R23
    sigma*I23

Rec_I24:
    I24 > R24
    sigma*I24

Rec_I31:
    I31 > R13
    sigma*I31

Rec_I32:
    I32 > R23
    sigma*I32

Rec_I34:
    I34 > R34
    sigma*I34

Rec_I41:
    I41 > R14
    sigma*I41

Rec_I42:
    I42 > R24
    sigma*I42

Rec_I43:
    I43 > R34
    sigma*I43
    
# Tertiary infections
# By DENV1
Inf_I231:
    R23 > I231
    R23*gamma*delta*beta*(I1+(phi*I21)+(phi*I31)+(phi*I41)+(phi*I231)+(phi*I241)+(phi*I341)+(phi*I2341))

Inf_I241:
    R24 > I241
    R24*gamma*delta*beta*(I1+(phi*I21)+(phi*I31)+(phi*I41)+(phi*I231)+(phi*I241)+(phi*I341)+(phi*I2341))

Inf_I341:
    R34 > I341
    R34*gamma*delta*beta*(I1+(phi*I21)+(phi*I31)+(phi*I41)+(phi*I231)+(phi*I241)+(phi*I341)+(phi*I2341))

# By DENV2
Inf_I132:
    R13 > I132
    R13*gamma*delta*beta*(I2+(phi*I12)+(phi*I32)+(phi*I42)+(phi*I132)+(phi*I142)+(phi*I342)+(phi*I1342))

Inf_I142:
    R14 > I142
    R14*gamma*delta*beta*(I2+(phi*I12)+(phi*I32)+(phi*I42)+(phi*I132)+(phi*I142)+(phi*I342)+(phi*I1342))

Inf_I342:
    R34 > I342
    R34*gamma*delta*beta*(I2+(phi*I12)+(phi*I32)+(phi*I42)+(phi*I132)+(phi*I142)+(phi*I342)+(phi*I1342))

# By DENV3
Inf_I123:
    R12 > I123
    R12*gamma*delta*beta*(I3+(phi*I13)+(phi*I23)+(phi*I43)+(phi*I123)+(phi*I143)+(phi*I243)+(phi*I1243))

Inf_I143:
    R14 > I143
    R14*gamma*delta*beta*(I3+(phi*I13)+(phi*I23)+(phi*I43)+(phi*I123)+(phi*I143)+(phi*I243)+(phi*I1243))

Inf_I243:
    R24 > I243
    R24*gamma*delta*beta*(I3+(phi*I13)+(phi*I23)+(phi*I43)+(phi*I123)+(phi*I143)+(phi*I243)+(phi*I1243))

#By DENV4
Inf_I124:
    R12 > I124
    R12*gamma*delta*beta*(I4+(phi*I14)+(phi*I24)+(phi*I34)+(phi*I124)+(phi*I134)+(phi*I234)+(phi*I1234))

Inf_I134:
    R13 > I134
    R13*gamma*delta*beta*(I4+(phi*I14)+(phi*I24)+(phi*I34)+(phi*I124)+(phi*I134)+(phi*I234)+(phi*I1234))

Inf_I234:
    R23 > I234
    R23*gamma*delta*beta*(I4+(phi*I14)+(phi*I24)+(phi*I34)+(phi*I124)+(phi*I134)+(phi*I234)+(phi*I1234))

# Tertiary Recoveries
# From DENV1
Rec_I231:
    I231 > R123
    sigma*I231

Rec_I241:
    I241 > R124
    sigma*I241

Rec_I341:
    I341 > R134
    sigma*I341
    
# From DENV2

Rec_I132:
    I132 > R123
    sigma*I132

Rec_I142:
    I142 > R124
    sigma*I142

Rec_I342:
    I342 > R234
    sigma*I342

# From DENV3
Rec_I123:
    I123 > R123
    sigma*I123

Rec_I143:
    I143 > R134
    sigma*I143

Rec_I243:
    I243 > R234
    sigma*I243

# From DENV4

Rec_I124:
    I124 > R124
    sigma*I124

Rec_I134:
    I134 > R134
    sigma*I134

Rec_I234:
    I234 > R234
    sigma*I234
    
# Quaternary Infections

Inf_I1234:
    R123 > I1234
    R123*gamma*delta*beta*(I4+(phi*I14)+(phi*I24)+(phi*I34)+(phi*I124)+(phi*I134)+(phi*I234)+(phi*I1234))

Inf_I1243:
    R124 > I1243
    R124*gamma*delta*beta*(I3+(phi*I13)+(phi*I23)+(phi*I43)+(phi*I123)+(phi*I143)+(phi*I243)+(phi*I1243))
    
Inf_I1342:
    R134 > I1342
    R134*gamma*delta*beta*(I2+(phi*I12)+(phi*I32)+(phi*I42)+(phi*I132)+(phi*I142)+(phi*I342)+(phi*I1342))
    
Inf_I2341:
    R234 > I2341
    R234*gamma*delta*beta*(I1+(phi*I21)+(phi*I31)+(phi*I41)+(phi*I231)+(phi*I241)+(phi*I341)+(phi*I2341))
  
# Quaternary Recoveries

Rec_I1234:
    I1234 > R1234
    sigma*I1234

Rec_I1243:
    I1243 > R1234
    sigma*I1243

Rec_I1342:
    I1342 > R1234
    sigma*I1342

Rec_I2341:
    I2341 > R1234
    sigma*I2341
    

#InitVar
S = 10
I1 = 50
I2 = 50
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
R12 = 0
R13 = 0
R14 = 0
R23 = 0
R24 = 0
R34 = 0
I231 = 0
I241 = 0
I341 = 0
I132 = 0
I142 = 0
I342 = 0
I123 = 0
I143 = 0
I243 = 0
I124 = 0
I134 = 0
I234 = 0
R123 = 0
R124 = 0
R134 = 0
R234 = 0
I1234 = 0
I1243 = 0
I1342 = 0
I2341 = 0
R1234 = 48000

#InitPar
mu = 1.0/(70*52) #70 years in weeks
beta = 400/52.0
phi = 1
sigma = 1/1.5
gamma = 1
delta = 0.2


#Event: denv4, _TIME_ > 40 , 0 {
#I4 = 50
#}
