
import testvis as tv

lam=3e8/(100e9)
#b1=tv.getbaselines('c40-1n.cfg',lam=lam)
#b2=tv.getbaselines('c40-2n.cfg',lam=lam)
#b3=tv.getbaselines('c40-3n.cfg',lam=lam)
#b4=tv.getbaselines('c40-4n.cfg',lam=lam)
#b5=tv.getbaselines('c40-5n.cfg',lam=lam)
#b6=tv.getbaselines('c40-6n.cfg',lam=lam)
#b7=tv.getbaselines('c40-7n2.cfg',lam=lam)
#b8=tv.getbaselines('c40-8n2.cfg',lam=lam)
#b9=tv.getbaselines('c40-9n2.cfg',lam=lam)
#baca=tv.getbaselines('aca.std10.cfg',lam=lam)

filelist=['aca.std10.cfg','c40-1n.cfg','c40-2n.cfg','c40-3n.cfg','c40-4n.cfg','c40-5n.cfg','c40-6n.cfg',
          'c40-7n2.cfg','c40-8n2.cfg','c40-9n2.cfg']

for ff in filelist:
    b=tv.getbaselines(ff,lam=lam)
    lasmin= 0.5/(b['qspecial'])[0]*206265.
    las5= 0.5/(b['qspecial'])[1]*206265.
    las10= 0.5/(b['qspecial'])[2]*206265.
    print ff,' ',lasmin,' ',las5,' ',las10

### RESULTS - 100GHz; 0.5 lambda/B_min , 0.5 lam/B_5 , 0.5 lam/B_10
#
#aca.std10.cfg   34.8   34.1   32.2
#c40-1n.cfg   20.6   14.6   11.2
#c40-2n.cfg   20.5   11.2   8.0
#c40-3n.cfg   20.5   6.9   5.1
#c40-4n.cfg   20.5   4.5   3.2
#c40-5n.cfg   18.5   3.0   1.9
#c40-6n.cfg   20.2   1.56   1.11
#c40-7n2.cfg   3.80   0.82   0.65
#c40-8n2.cfg   1.84   0.62   0.45
#c40-9n2.cfg   1.14   0.44   0.30
#
