program tzeta

!   This tests the special function routines:
!     mpberne, besselj, gamma, incgamma, zeta and zetaem.
!   David H Bailey   22 May 2015

use mpmodule
implicit none
integer nb, ndp, nwds
parameter (nb = 1024, ndp = 1000, nwds = int (ndp / mpdpw + 2))
double precision d1, d2, second, tm0, tm1, tm2
type (mp_real) berne(nb), t1, t2, t3
external second

tm2 = second ()
write (6, 1) nb
1 format ('Generate even Bernoulli coefficient array; nb =',i8/ &
  'berne(1), berne(nb) =')
tm0 = second ()
call mpberne (nb, berne, nwds)
tm1 = second ()
tm2 = tm1 - tm0
call mpwrite (6, ndp + 20, ndp, berne(1), berne(nb))
write (6, 2) tm1 - tm0
2 format ('CPU time =',f10.4)

write (6, 11)
11 format ('besselj (2, 12) =')
t1 = mpreal (2.d0, nwds)
t2 = mpreal (12.d0, nwds)
tm0 = second ()
t3 = besselj (t1, t2)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
call mpwrite (6, ndp + 20, ndp, t3)
write (6, 2) tm1 - tm0

write (6, 12)
12 format ('besselj (2.5, 12.125) =')
t1 = mpreal (2.5d0, nwds)
t2 = mpreal (12.125d0, nwds)
tm0 = second ()
t3 = besselj (t1, t2)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
call mpwrite (6, ndp + 20, ndp, t3)
write (6, 2) tm1 - tm0

write (6, 21)
21 format ('gamma (7.25) =')
t1 = mpreal (7.25d0, nwds)
tm0 = second ()
t3 = gamma (t1)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
call mpwrite (6, ndp + 20, ndp, t3)
write (6, 2) tm1 - tm0

write (6, 22)
22 format ('gamma (pi) =')
t1 = mppi (nwds)
tm0 = second ()
t3 = gamma (t1)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
call mpwrite (6, ndp + 20, ndp, t3)
write (6, 2) tm1 - tm0

write (6, 31) 
31 format ('incgamma (2, 12) =')
t1 = mpreal (2.d0, nwds)
t2 = mpreal (12.d0, nwds)
tm0 = second ()
t3 = incgamma (t1, t2)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
call mpwrite (6, ndp + 20, ndp, t3)
write (6, 2) tm1 - tm0

write (6, 32)
32 format ('incgamma (2.5, 12.125) =')
t1 = mpreal (2.5d0, nwds)
t2 = mpreal (12.125d0, nwds)
tm0 = second ()
t3 = incgamma (t1, t2)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
call mpwrite (6, ndp + 20, ndp, t3)
write (6, 2) tm1 - tm0

write (6, 41)
41 format ('zeta(3) =')
t1 = mpreal (3.d0, nwds)
tm0 = second ()
t2 = zeta (t1)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
call mpwrite (6, ndp + 20, ndp, t2)
write (6, 2) tm1 - tm0

write (6, 42)
42 format ('zeta(pi) =')
t1 = mppi (nwds)
tm0 = second ()
t2 = zeta (t1)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
call mpwrite (6, ndp + 20, ndp, t2)
write (6, 2) tm1 - tm0

write (6, 43)
43 format ('zeta(-4.125) =')
t1 = mpreal (-4.125d0, nwds)
tm0 = second ()
t2 = zeta (t1)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
call mpwrite (6, ndp + 20, ndp, t2)
write (6, 2) tm1 - tm0

write (6, 44)
44 format ('zeta(-13) =')
t1 = mpreal (-13.d0, nwds)
tm0 = second ()
t2 = zeta (t1)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
call mpwrite (6, ndp + 20, ndp, t2)
write (6, 2) tm1 - tm0

write (6, 45)
45 format ('zeta(1033) =')
t1 = mpreal (1033.d0, nwds)
tm0 = second ()
t2 = zeta (t1)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
call mpwrite (6, ndp + 20, ndp, t2)
write (6, 2) tm1 - tm0

write (6, 51)
51 format ('zetaem(3) =')
t1 = mpreal (3.d0, nwds)
tm0 = second ()
t2 = zetaem (nb, berne, t1)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
call mpwrite (6, ndp + 20, ndp, t2)
write (6, 2) tm1 - tm0

write (6, 52)
52 format ('zetaem(pi) =')
t1 = mppi (nwds)
tm0 = second ()
t2 = zetaem (nb, berne, t1)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
call mpwrite (6, ndp + 20, ndp, t2)
write (6, 2) tm1 - tm0

write (6, 53)
53 format ('zetaem(-4.125) =')
t1 = mpreal (-4.125d0, nwds)
tm0 = second ()
t2 = zetaem (nb, berne, t1)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
call mpwrite (6, ndp + 20, ndp, t2)
write (6, 2) tm1 - tm0

write (6, 54)
54 format ('zetaem(-13) =')
t1 = mpreal (-13.d0, nwds)
tm0 = second ()
t2 = zetaem (nb, berne, t1)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
call mpwrite (6, ndp + 20, ndp, t2)
write (6, 2) tm1 - tm0

write (6, 55)
55 format ('zetaem(1033) =')
t1 = mpreal (1033.d0, nwds)
tm0 = second ()
t2 = zetaem (nb, berne, t1)
tm1 = second ()
tm2 = tm2 + (tm1 - tm0)
call mpwrite (6, ndp + 20, ndp, t2)
write (6, 2) tm1 - tm0

write (6, 999) tm2
999 format ('Total CPU time =',f10.4)

stop
end


