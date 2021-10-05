C*********************************************************************
C*** FORTRAN 77 subroutines that call C functions in MathRC.c
C*********************************************************************

       subroutine getRandgenStateRF77()
       call getrandgenstaterc()
       end

       subroutine putRandgenStateRF77()
       call putrandgenstaterc()
       end

C*********************************************************************

       subroutine dunifRF77(x, a, b, logd, ans)
       double precision x, a, b, ans
       integer logd
       call dunifrc(x, a, b, logd, ans)
       end

       subroutine punifRF77(q, a, b, lowertail, logp, ans)
       double precision q, a, b, ans
       integer lowertail, logp
       call punifrc(q, a, b, lowertail, logp, ans)
       end

       subroutine qunifRF77(p, a, b, lowertail, logp, ans)
       double precision p, a, b, ans
       integer lowertail, logp
       call qunifrc(p, a, b, lowertail, logp, ans)
       end

       subroutine runifRF77(a, b, ans)
       double precision a, b, ans
       call runifrc(a, b, ans)
       end

C*********************************************************************

       subroutine dnormRF77(x, mu, sigma, logd, ans)
       double precision x, mu, sigma, ans
       integer logd
       call dnormrc(x, mu, sigma, logd, ans)
       end

       subroutine pnormRF77(q, mu, sigma, lowertail, logp, ans)
       double precision q, mu, sigma, ans
       integer lowertail, logp
       call pnormrc(q, mu, sigma, lowertail, logp, ans)
       end

       subroutine qnormRF77(p, mu, sigma, lowertail, logp, ans)
       double precision p, mu, sigma, ans
       integer lowertail, logp
       call qnormrc(p, mu, sigma, lowertail, logp, ans)
       end

       subroutine rnormRF77(mu, sigma, ans)
       double precision mu, sigma, ans
       call rnormrc(mu, sigma, ans)
       end

C*********************************************************************

       subroutine dchisqRF77(x, df, logd, ans)
       double precision x, df, ans
       integer logd
       call dchisqrc(x, df, logd, ans)
       end

       subroutine pchisqRF77(q, df, lowertail, logp, ans)
       double precision q, df, ans
       integer lowertail, logp
       call pchisqrc(q, df, lowertail, logp, ans)
       end

       subroutine qchisqRF77(p, df, lowertail, logp, ans)
       double precision p, df, ans
       integer lowertail, logp
       call qchisqrc(p, df, lowertail, logp, ans)
       end

       subroutine rchisqRF77(df, ans)
       double precision df, ans
       call rchisqrc(df, ans)
       end

C*********************************************************************

       subroutine dtRF77(x, df, logd, ans)
       double precision x, df, ans
       integer logd
       call dtrc(x, df, logd, ans)
       end

       subroutine ptRF77(q, df, lowertail, logp, ans)
       double precision q, df, ans
       integer lowertail, logp
       call ptrc(q, df, lowertail, logp, ans)
       end

       subroutine qtRF77(p, df, lowertail, logp, ans)
       double precision p, df, ans
       integer lowertail, logp
       call qtrc(p, df, lowertail, logp, ans)
       end

       subroutine rtRF77(df, ans)
       double precision df, ans
       call rtrc(df, ans)
       end

C*********************************************************************

       subroutine dbetaRF77(x, shape1, shape2, logd, ans)
       double precision x, shape1, shape2, ans
       integer logd
       call dbetarc(x, shape1, shape2, logd, ans)
       end

       subroutine pbetaRF77(q, shape1, shape2, lowertail, logp, ans)
       double precision q, shape1, shape2, ans
       integer lowertail, logp
       call pbetarc(q, shape1, shape2, lowertail, logp, ans)
       end

       subroutine qbetaRF77(p, shape1, shape2, lowertail, logp, ans)
       double precision p, shape1, shape2, ans
       integer lowertail, logp
       call qbetarc(p, shape1, shape2, lowertail, logp, ans)
       end

       subroutine rbetaRF77(shape1, shape2, ans)
       double precision shape1, shape2, ans
       call rbetarc(shape1, shape2, ans)
       end

C*********************************************************************

       subroutine dgammaRF77(x, shape, scale, logd, ans)
       double precision x, shape, scale, ans
       integer logd
       call dgammarc(x, shape, scale, logd, ans)
       end

       subroutine pgammaRF77(q, shape, scale, lowertail, logp, ans)
       double precision q, shape, scale, ans
       integer lowertail, logp
       call pgammarc(q, shape, scale, lowertail, logp, ans)
       end

       subroutine qgammaRF77(p, shape, scale, lowertail, logp, ans)
       double precision p, shape, scale, ans
       integer lowertail, logp
       call qgammarc(p, shape, scale, lowertail, logp, ans)
       end

       subroutine rgammaRF77(shape, scale, ans)
       double precision shape, scale, ans
       call rgammarc(shape, scale, ans)
       end

C*********************************************************************

       subroutine dbinomRF77(x, n, prob, logd, ans)
       double precision x, n, prob, ans
       integer logd
       call dbinomrc(x, n, prob, logd, ans)
       end

       subroutine dbinomrawRF77(x, n, probp, probq, logd, ans)
       double precision x, n, probp, probq, ans
       integer logd
       call dbinomrawrc(x, n, probp, probq, logd, ans)
       end

       subroutine pbinomRF77(q, n, prob, lowertail, logp, ans)
       double precision q, n, prob, ans
       integer lowertail, logp
       call pbinomrc(q, n, prob, lowertail, logp, ans)
       end

       subroutine qbinomRF77(p, n, prob, lowertail, logp, ans)
       double precision p, n, prob, ans
       integer lowertail, logp
       call qbinomrc(p, n, prob, lowertail, logp, ans)
       end

       subroutine rbinomRF77(n, prob, ans)
       double precision n, prob, ans
       call rbinomrc(n, prob, ans)
       end

C*********************************************************************

       subroutine dpoisRF77(x, lambda, logd, ans)
       double precision x, lambda, ans
       integer logd
       call dpoisrc(x, lambda, logd, ans)
       end

       subroutine dpoisrawRF77(x, lambda, logd, ans)
       double precision x, lambda, ans
       integer logd
       call dpoisrawrc(x, lambda, logd, ans)
       end

       subroutine ppoisRF77(q, lambda, lowertail, logp, ans)
       double precision q, lambda, ans
       integer lowertail, logp
       call ppoisrc(q, lambda, lowertail, logp, ans)
       end

       subroutine qpoisRF77(p, lambda, lowertail, logp, ans)
       double precision p, lambda, ans
       integer lowertail, logp
       call qpoisrc(p, lambda, lowertail, logp, ans)
       end

       subroutine rpoisRF77(lambda, ans)
       double precision lambda, ans
       call rpoisrc(lambda, ans)
       end

C*********************************************************************

       subroutine dnbinomRF77(x, n, prob, logd, ans)
       double precision x, n, prob, ans
       integer logd
       call dnbinomrc(x, n, prob, logd, ans)
       end

       subroutine pnbinomRF77(q, n, prob, lowertail, logp, ans)
       double precision q, n, prob, ans
       integer lowertail, logp
       call pnbinomrc(q, n, prob, lowertail, logp, ans)
       end

       subroutine qnbinomRF77(p, n, prob, lowertail, logp, ans)
       double precision p, n, prob, ans
       integer lowertail, logp
       call qnbinomrc(p, n, prob, lowertail, logp, ans)
       end

       subroutine rnbinomRF77(n, prob, ans)
       double precision n, prob, ans
       call rnbinomrc(n, prob, ans)
       end

C*********************************************************************

       subroutine lgammaRF77(x, ans)
       double precision x, ans
       call lgammarc(x, ans)
       end

       subroutine digammaRF77(x, ans)
       double precision x, ans
       call digammarc(x, ans)
       end

       subroutine trigammaRF77(x, ans)
       double precision x, ans
       call trigammarc(x, ans)
       end

C*********************************************************************

       subroutine log1pRF77(x, ans)
       double precision x, ans
       call log1prc(x, ans)
       end

       subroutine log1pmxRF77(x, ans)
       double precision x, ans
       call log1pmxrc(x, ans)
       end

       subroutine log1pexpRF77(x, ans)
       double precision x, ans
       call log1pexprc(x, ans)
       end

       subroutine expm1RF77(x, ans)
       double precision x, ans
       call expm1rc(x, ans)
       end

       subroutine lgamma1pRF77(x, ans)
       double precision x, ans
       call lgamma1prc(x, ans)
       end

C*********************************************************************
