!----------------------------------------------------------------------
!----------------------------------------------------------------------
!   AUTO file for qgs model
!----------------------------------------------------------------------
!----------------------------------------------------------------------

SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP)
	!--------- ----

	! Evaluates the algebraic equations or ODE right hand side

	! Input arguments :
	!      NDIM   :   Dimension of the algebraic or ODE system 
	!      U      :   State variables
	!      ICP    :   Array indicating the free parameter(s)
	!      PAR    :   Equation parameters

	! Values to be returned :
	!      F      :   Equation or ODE right hand side values
  
	! Normally unused Jacobian arguments : IJAC, DFDU, DFDP (see manual)

	IMPLICIT NONE
	INTEGER, INTENT(IN) :: NDIM, IJAC, ICP(*)
	DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
	DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
	DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM),DFDP(NDIM,*)

	DOUBLE PRECISION thetas_1
	DOUBLE PRECISION hk_2
	DOUBLE PRECISION k_d
	DOUBLE PRECISION k_p
	DOUBLE PRECISION sigma

	thetas_1 = PAR(1)
	hk_2 = PAR(2)
	k_d = PAR(3)
	k_p = PAR(4)
	sigma = PAR(5)


	F(1) =	 -0.5*U(1)*k_d + 0.5*U(11)*k_d - 0.780274140669492*U(13)*hk_2 + 0.780274140669492&
		&*U(3)*hk_2

	F(2) =	 -0.980418808722262*U(1)*U(3) - 0.980418808722262*U(11)*U(13) + 0.5*U(12)*k_d - 1&
		&.56867009395562*U(14)*U(16) - 1.50055762081784*U(15)*U(18) + 1.50055762081784*U(&
		&16)*U(17) - 0.5*U(2)*k_d + 0.101317695204044*U(3) - 1.56867009395562*U(4)*U(6) -&
		& 1.50055762081784*U(5)*U(8) + 1.50055762081784*U(6)*U(7)

	F(3) =	 0.980418808722261*U(1)*U(2) - 0.290064736308361*U(1)*hk_2 + 0.980418808722261*U(&
		&11)*U(12) + 0.290064736308361*U(11)*hk_2 + 0.5*U(13)*k_d + 1.56867009395562*U(14&
		&)*U(15) + 1.50055762081784*U(15)*U(17) + 1.50055762081784*U(16)*U(18) - 0.101317&
		&695204044*U(2) - 0.5*U(3)*k_d + 1.56867009395562*U(4)*U(5) + 1.50055762081784*U(&
		&5)*U(7) + 1.50055762081784*U(6)*U(8)

	F(4) =	 3.74531587521356*U(10)*U(7) + 1.87265793760678*U(12)*U(16) - 1.87265793760678*U(&
		&13)*U(15) + 0.5*U(14)*k_d - 0.312109656267797*U(16)*hk_2 + 3.74531587521356*U(17&
		&)*U(20) - 3.74531587521356*U(18)*U(19) + 1.87265793760678*U(2)*U(6) - 1.87265793&
		&760678*U(3)*U(5) - 0.5*U(4)*k_d + 0.312109656267797*U(6)*hk_2 - 3.74531587521356&
		&*U(8)*U(9)

	F(5) =	 -1.02902937637678*U(1)*U(6) - 1.02902937637678*U(11)*U(16) + 1.73752196836555*U(&
		&12)*U(18) + 0.574852231579352*U(13)*U(14) - 1.73752196836555*U(13)*U(17) + 0.5*U&
		&(15)*k_d - 0.171353251318102*U(18)*hk_2 + 1.73752196836555*U(2)*U(8) + 0.5748522&
		&31579352*U(3)*U(4) - 1.73752196836555*U(3)*U(7) - 0.5*U(5)*k_d + 0.0478988752370&
		&612*U(6) + 0.171353251318102*U(8)*hk_2

	F(6) =	 1.02902937637678*U(1)*U(5) + 1.02902937637678*U(11)*U(15) - 0.574852231579352*U(&
		&12)*U(14) - 1.73752196836555*U(12)*U(17) - 1.73752196836555*U(13)*U(18) + 0.2194&
		&09248694409*U(14)*hk_2 + 0.5*U(16)*k_d + 0.171353251318102*U(17)*hk_2 - 0.574852&
		&231579352*U(2)*U(4) - 1.73752196836555*U(2)*U(7) - 1.73752196836555*U(3)*U(8) - &
		&0.219409248694409*U(4)*hk_2 - 0.0478988752370612*U(5) - 0.5*U(6)*k_d - 0.1713532&
		&51318102*U(7)*hk_2

	F(7) =	 -2.71889339738442*U(1)*U(8) - 4.35022943581506*U(10)*U(4) - 2.71889339738442*U(1&
		&1)*U(18) + 0.753865979381444*U(12)*U(16) + 0.753865979381444*U(13)*U(15) - 4.350&
		&22943581506*U(14)*U(20) - 0.125644329896907*U(16)*hk_2 + 0.5*U(17)*k_d + 0.75386&
		&5979381444*U(2)*U(6) + 0.753865979381444*U(3)*U(5) + 0.125644329896907*U(6)*hk_2&
		& - 0.5*U(7)*k_d + 0.0702434536337316*U(8)

	F(8) =	 2.71889339738442*U(1)*U(7) + 2.71889339738442*U(11)*U(17) - 0.753865979381443*U(&
		&12)*U(15) + 0.753865979381443*U(13)*U(16) + 4.35022943581507*U(14)*U(19) + 0.125&
		&644329896907*U(15)*hk_2 + 0.5*U(18)*k_d - 0.753865979381443*U(2)*U(5) + 0.753865&
		&979381443*U(3)*U(6) + 4.35022943581507*U(4)*U(9) - 0.125644329896907*U(5)*hk_2 -&
		& 0.0702434536337315*U(7) - 0.5*U(8)*k_d

	F(9) =	 -2.26482546109568*U(1)*U(10) + 0.050658847602022*U(10) - 2.26482546109568*U(11)*&
		&U(20) - 1.7450294536311*U(14)*U(18) + 0.5*U(19)*k_d - 1.7450294536311*U(4)*U(8) &
		&- 0.5*U(9)*k_d

	F(10) =	 2.26482546109569*U(1)*U(9) - 0.5*U(10)*k_d + 2.26482546109569*U(11)*U(19) + 4.12&
		&72231398691e-17*U(12)**2 - 4.1272231398691e-17*U(13)**2 + 1.7450294536311*U(14)*&
		&U(17) + 4.1272231398691e-17*U(2)**2 + 0.5*U(20)*k_d - 4.1272231398691e-17*U(3)**&
		&2 + 1.7450294536311*U(4)*U(7) - 0.050658847602022*U(9)

	F(11) =	 0.25*U(1)*k_d*sigma/(0.5*sigma + 1.0) - 2.49687725014237*U(10)*U(19)/(0.5*sigma &
		&+ 1.0) + U(11)*(-0.5*sigma*(0.5*k_d + 2.0*k_p)/(0.5*sigma + 1.0) - 0.045/(0.5*si&
		&gma + 1.0)) - 1.56054828133898*U(12)*U(3)/(0.5*sigma + 1.0) + 1.56054828133898*U&
		&(13)*U(2)/(0.5*sigma + 1.0) + 0.390137070334746*U(13)*hk_2*sigma/(0.5*sigma + 1.&
		&0) - 1.24843862507119*U(15)*U(6)/(0.5*sigma + 1.0) + 1.24843862507119*U(16)*U(5)&
		&/(0.5*sigma + 1.0) - 3.12109656267797*U(17)*U(8)/(0.5*sigma + 1.0) + 3.121096562&
		&67797*U(18)*U(7)/(0.5*sigma + 1.0) + 2.49687725014237*U(20)*U(9)/(0.5*sigma + 1.&
		&0) - 0.390137070334746*U(3)*hk_2*sigma/(0.5*sigma + 1.0) + 0.045*thetas_1/(0.5*s&
		&igma + 1.0)

	F(12) =	 U(1)*U(13)*(-1.31866329773144*sigma/(1.345*sigma + 1.0) - 1.56054828133898/(1.34&
		&5*sigma + 1.0)) + U(11)*U(3)*(-1.31866329773144*sigma/(1.345*sigma + 1.0) + 1.56&
		&054828133898/(1.345*sigma + 1.0)) + U(12)*(-1.345*sigma*(0.5*k_d + 2.0*k_p)/(1.3&
		&45*sigma + 1.0) - 0.045/(1.345*sigma + 1.0)) + 0.136272300049439*U(13)*sigma/(1.&
		&345*sigma + 1.0) + U(14)*U(6)*(-2.10986127637031*sigma/(1.345*sigma + 1.0) + 2.4&
		&9687725014237/(1.345*sigma + 1.0)) + U(15)*U(8)*(-2.01825*sigma/(1.345*sigma + 1&
		&.0) + 1.95/(1.345*sigma + 1.0)) + U(16)*U(4)*(-2.10986127637031*sigma/(1.345*sig&
		&ma + 1.0) - 2.49687725014237/(1.345*sigma + 1.0)) + U(16)*U(7)*(2.01825*sigma/(1&
		&.345*sigma + 1.0) - 1.95/(1.345*sigma + 1.0)) + U(17)*U(6)*(2.01825*sigma/(1.345&
		&*sigma + 1.0) + 1.95/(1.345*sigma + 1.0)) + U(18)*U(5)*(-2.01825*sigma/(1.345*si&
		&gma + 1.0) - 1.95/(1.345*sigma + 1.0)) + 0.6725*U(2)*k_d*sigma/(1.345*sigma + 1.&
		&0)

	F(13) =	 U(1)*U(12)*(1.31866329773144*sigma/(1.345*sigma + 1.0) + 1.56054828133898/(1.345&
		&*sigma + 1.0)) + 0.390137070334746*U(1)*hk_2*sigma/(1.345*sigma + 1.0) + U(11)*U&
		&(2)*(1.31866329773144*sigma/(1.345*sigma + 1.0) - 1.56054828133898/(1.345*sigma &
		&+ 1.0)) - 0.390137070334746*U(11)*hk_2*sigma/(1.345*sigma + 1.0) - 0.13627230004&
		&9439*U(12)*sigma/(1.345*sigma + 1.0) + U(13)*(-1.345*sigma*(0.5*k_d + 2.0*k_p)/(&
		&1.345*sigma + 1.0) - 0.045/(1.345*sigma + 1.0)) + U(14)*U(5)*(2.10986127637031*s&
		&igma/(1.345*sigma + 1.0) - 2.49687725014237/(1.345*sigma + 1.0)) + U(15)*U(4)*(2&
		&.10986127637031*sigma/(1.345*sigma + 1.0) + 2.49687725014237/(1.345*sigma + 1.0)&
		&) + U(15)*U(7)*(2.01825*sigma/(1.345*sigma + 1.0) - 1.95/(1.345*sigma + 1.0)) + &
		&U(16)*U(8)*(2.01825*sigma/(1.345*sigma + 1.0) - 1.95/(1.345*sigma + 1.0)) + U(17&
		&)*U(5)*(2.01825*sigma/(1.345*sigma + 1.0) + 1.95/(1.345*sigma + 1.0)) + U(18)*U(&
		&6)*(2.01825*sigma/(1.345*sigma + 1.0) + 1.95/(1.345*sigma + 1.0)) + 0.6725*U(3)*&
		&k_d*sigma/(1.345*sigma + 1.0)

	F(14) =	 U(10)*U(17)*(7.49063175042712*sigma/(2.0*sigma + 1.0) - 4.99375450028475/(2.0*si&
		&gma + 1.0)) + U(12)*U(6)*(3.74531587521356*sigma/(2.0*sigma + 1.0) - 2.496877250&
		&14237/(2.0*sigma + 1.0)) + U(13)*U(5)*(-3.74531587521356*sigma/(2.0*sigma + 1.0)&
		& + 2.49687725014237/(2.0*sigma + 1.0)) + U(14)*(-2.0*sigma*(0.5*k_d + 2.0*k_p)/(&
		&2.0*sigma + 1.0) - 0.045/(2.0*sigma + 1.0)) + U(15)*U(3)*(-3.74531587521356*sigm&
		&a/(2.0*sigma + 1.0) - 2.49687725014237/(2.0*sigma + 1.0)) + U(16)*U(2)*(3.745315&
		&87521356*sigma/(2.0*sigma + 1.0) + 2.49687725014237/(2.0*sigma + 1.0)) + 0.62421&
		&9312535594*U(16)*hk_2*sigma/(2.0*sigma + 1.0) + U(18)*U(9)*(-7.49063175042712*si&
		&gma/(2.0*sigma + 1.0) + 4.99375450028475/(2.0*sigma + 1.0)) + U(19)*U(8)*(-7.490&
		&63175042712*sigma/(2.0*sigma + 1.0) - 4.99375450028475/(2.0*sigma + 1.0)) + U(20&
		&)*U(7)*(7.49063175042712*sigma/(2.0*sigma + 1.0) + 4.99375450028475/(2.0*sigma +&
		& 1.0)) + 1.0*U(4)*k_d*sigma/(2.0*sigma + 1.0) - 0.624219312535594*U(6)*hk_2*sigm&
		&a/(2.0*sigma + 1.0)

	F(15) =	 U(1)*U(16)*(-2.92758857579193*sigma/(2.845*sigma + 1.0) - 1.24843862507119/(2.84&
		&5*sigma + 1.0)) + U(11)*U(6)*(-2.92758857579193*sigma/(2.845*sigma + 1.0) + 1.24&
		&843862507119/(2.845*sigma + 1.0)) + U(12)*U(8)*(4.94325*sigma/(2.845*sigma + 1.0&
		&) - 1.95/(2.845*sigma + 1.0)) + U(13)*U(4)*(1.63545459884326*sigma/(2.845*sigma &
		&+ 1.0) - 2.49687725014237/(2.845*sigma + 1.0)) + U(13)*U(7)*(-4.94325*sigma/(2.8&
		&45*sigma + 1.0) + 1.95/(2.845*sigma + 1.0)) + U(14)*U(3)*(1.63545459884326*sigma&
		&/(2.845*sigma + 1.0) + 2.49687725014237/(2.845*sigma + 1.0)) + U(15)*(-2.845*sig&
		&ma*(0.5*k_d + 2.0*k_p)/(2.845*sigma + 1.0) - 0.045/(2.845*sigma + 1.0)) + 0.1362&
		&72300049439*U(16)*sigma/(2.845*sigma + 1.0) + U(17)*U(3)*(-4.94325*sigma/(2.845*&
		&sigma + 1.0) - 1.95/(2.845*sigma + 1.0)) + U(18)*U(2)*(4.94325*sigma/(2.845*sigm&
		&a + 1.0) + 1.95/(2.845*sigma + 1.0)) + 0.4875*U(18)*hk_2*sigma/(2.845*sigma + 1.&
		&0) + 1.4225*U(5)*k_d*sigma/(2.845*sigma + 1.0) - 0.4875*U(8)*hk_2*sigma/(2.845*s&
		&igma + 1.0)

	F(16) =	 U(1)*U(15)*(2.92758857579193*sigma/(2.845*sigma + 1.0) + 1.24843862507119/(2.845&
		&*sigma + 1.0)) + U(11)*U(5)*(2.92758857579193*sigma/(2.845*sigma + 1.0) - 1.2484&
		&3862507119/(2.845*sigma + 1.0)) + U(12)*U(4)*(-1.63545459884326*sigma/(2.845*sig&
		&ma + 1.0) + 2.49687725014237/(2.845*sigma + 1.0)) + U(12)*U(7)*(-4.94325*sigma/(&
		&2.845*sigma + 1.0) + 1.95/(2.845*sigma + 1.0)) + U(13)*U(8)*(-4.94325*sigma/(2.8&
		&45*sigma + 1.0) + 1.95/(2.845*sigma + 1.0)) + U(14)*U(2)*(-1.63545459884326*sigm&
		&a/(2.845*sigma + 1.0) - 2.49687725014237/(2.845*sigma + 1.0)) - 0.62421931253559&
		&4*U(14)*hk_2*sigma/(2.845*sigma + 1.0) - 0.136272300049439*U(15)*sigma/(2.845*si&
		&gma + 1.0) + U(16)*(-2.845*sigma*(0.5*k_d + 2.0*k_p)/(2.845*sigma + 1.0) - 0.045&
		&/(2.845*sigma + 1.0)) + U(17)*U(2)*(-4.94325*sigma/(2.845*sigma + 1.0) - 1.95/(2&
		&.845*sigma + 1.0)) - 0.4875*U(17)*hk_2*sigma/(2.845*sigma + 1.0) + U(18)*U(3)*(-&
		&4.94325*sigma/(2.845*sigma + 1.0) - 1.95/(2.845*sigma + 1.0)) + 0.62421931253559&
		&4*U(4)*hk_2*sigma/(2.845*sigma + 1.0) + 1.4225*U(6)*k_d*sigma/(2.845*sigma + 1.0&
		&) + 0.4875*U(7)*hk_2*sigma/(2.845*sigma + 1.0)

	F(17) =	 U(1)*U(18)*(-10.5493063818515*sigma/(3.88*sigma + 1.0) - 3.12109656267797/(3.88*&
		&sigma + 1.0)) + U(10)*U(14)*(-16.8788902109624*sigma/(3.88*sigma + 1.0) + 4.9937&
		&5450028475/(3.88*sigma + 1.0)) + U(11)*U(8)*(-10.5493063818515*sigma/(3.88*sigma&
		& + 1.0) + 3.12109656267797/(3.88*sigma + 1.0)) + U(12)*U(6)*(2.925*sigma/(3.88*s&
		&igma + 1.0) - 1.95/(3.88*sigma + 1.0)) + U(13)*U(5)*(2.925*sigma/(3.88*sigma + 1&
		&.0) - 1.95/(3.88*sigma + 1.0)) + U(15)*U(3)*(2.925*sigma/(3.88*sigma + 1.0) + 1.&
		&95/(3.88*sigma + 1.0)) + U(16)*U(2)*(2.925*sigma/(3.88*sigma + 1.0) + 1.95/(3.88&
		&*sigma + 1.0)) + 0.4875*U(16)*hk_2*sigma/(3.88*sigma + 1.0) + U(17)*(-3.88*sigma&
		&*(0.5*k_d + 2.0*k_p)/(3.88*sigma + 1.0) - 0.045/(3.88*sigma + 1.0)) + 0.27254460&
		&0098878*U(18)*sigma/(3.88*sigma + 1.0) + U(20)*U(4)*(-16.8788902109624*sigma/(3.&
		&88*sigma + 1.0) - 4.99375450028475/(3.88*sigma + 1.0)) - 0.4875*U(6)*hk_2*sigma/&
		&(3.88*sigma + 1.0) + 1.94*U(7)*k_d*sigma/(3.88*sigma + 1.0)

	F(18) =	 U(1)*U(17)*(10.5493063818515*sigma/(3.88*sigma + 1.0) + 3.12109656267797/(3.88*s&
		&igma + 1.0)) + U(11)*U(7)*(10.5493063818515*sigma/(3.88*sigma + 1.0) - 3.1210965&
		&6267797/(3.88*sigma + 1.0)) + U(12)*U(5)*(-2.925*sigma/(3.88*sigma + 1.0) + 1.95&
		&/(3.88*sigma + 1.0)) + U(13)*U(6)*(2.925*sigma/(3.88*sigma + 1.0) - 1.95/(3.88*s&
		&igma + 1.0)) + U(14)*U(9)*(16.8788902109624*sigma/(3.88*sigma + 1.0) - 4.9937545&
		&0028475/(3.88*sigma + 1.0)) + U(15)*U(2)*(-2.925*sigma/(3.88*sigma + 1.0) - 1.95&
		&/(3.88*sigma + 1.0)) - 0.4875*U(15)*hk_2*sigma/(3.88*sigma + 1.0) + U(16)*U(3)*(&
		&2.925*sigma/(3.88*sigma + 1.0) + 1.95/(3.88*sigma + 1.0)) - 0.272544600098878*U(&
		&17)*sigma/(3.88*sigma + 1.0) + U(18)*(-3.88*sigma*(0.5*k_d + 2.0*k_p)/(3.88*sigm&
		&a + 1.0) - 0.045/(3.88*sigma + 1.0)) + U(19)*U(4)*(16.8788902109624*sigma/(3.88*&
		&sigma + 1.0) + 4.99375450028475/(3.88*sigma + 1.0)) + 0.4875*U(5)*hk_2*sigma/(3.&
		&88*sigma + 1.0) + 1.94*U(8)*k_d*sigma/(3.88*sigma + 1.0)

	F(19) =	 U(1)*U(20)*(-12.1847609806948*sigma/(5.38*sigma + 1.0) - 2.49687725014237/(5.38*&
		&sigma + 1.0)) + U(10)*U(11)*(-12.1847609806948*sigma/(5.38*sigma + 1.0) + 2.4968&
		&7725014237/(5.38*sigma + 1.0)) + U(14)*U(8)*(-9.38825846053533*sigma/(5.38*sigma&
		& + 1.0) + 4.99375450028475/(5.38*sigma + 1.0)) + U(18)*U(4)*(-9.38825846053533*s&
		&igma/(5.38*sigma + 1.0) - 4.99375450028475/(5.38*sigma + 1.0)) + U(19)*(-5.38*si&
		&gma*(0.5*k_d + 2.0*k_p)/(5.38*sigma + 1.0) - 0.045/(5.38*sigma + 1.0)) + 0.27254&
		&4600098878*U(20)*sigma/(5.38*sigma + 1.0) + 2.69*U(9)*k_d*sigma/(5.38*sigma + 1.&
		&0)

	F(20) =	 U(1)*U(19)*(12.1847609806948*sigma/(5.38*sigma + 1.0) + 2.49687725014237/(5.38*s&
		&igma + 1.0)) + 2.69*U(10)*k_d*sigma/(5.38*sigma + 1.0) + U(11)*U(9)*(12.18476098&
		&06948*sigma/(5.38*sigma + 1.0) - 2.49687725014237/(5.38*sigma + 1.0)) + 4.440892&
		&09849915e-16*U(12)*U(2)*sigma/(5.38*sigma + 1.0) - 4.44089209849915e-16*U(13)*U(&
		&3)*sigma/(5.38*sigma + 1.0) + U(14)*U(7)*(9.38825846053533*sigma/(5.38*sigma + 1&
		&.0) - 4.99375450028475/(5.38*sigma + 1.0)) + U(17)*U(4)*(9.38825846053533*sigma/&
		&(5.38*sigma + 1.0) + 4.99375450028475/(5.38*sigma + 1.0)) - 0.272544600098878*U(&
		&19)*sigma/(5.38*sigma + 1.0) + U(20)*(-5.38*sigma*(0.5*k_d + 2.0*k_p)/(5.38*sigm&
		&a + 1.0) - 0.045/(5.38*sigma + 1.0))


END SUBROUTINE FUNC

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

SUBROUTINE STPNT(NDIM,U,PAR,T)
	!--------- -----
  
	! Input arguments :
	!      NDIM   :   Dimension of the algebraic or ODE system 

	! Values to be returned :
	!      U      :   A starting solution vector
	!      PAR    :   The corresponding equation-parameter values

	! Note : For time- or space-dependent solutions this subroutine has
	!        the scalar input parameter T contains the varying time or space
	!        variable value.

	IMPLICIT NONE
	INTEGER, INTENT(IN) :: NDIM
	DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
	DOUBLE PRECISION, INTENT(IN) :: T
	DOUBLE PRECISION :: X(NDIM+1)
	INTEGER :: i,is

	! Initialize the equation parameters

! INITIALISE PARAMETERS

	! Initialize the solution

! INITIALISE SOLUTION

	! Initialization from a solution file (selection with PAR36)
	! open(unit=15,file='',status='old')
	! is=int(PAR(36))
	! if (is.gt.0) print*, 'Loading from solution :',is
	! DO i=1,is
	!    read(15,*) X
	! ENDDO
	! close(15)
	! U=X(2:NDIM+1)

END SUBROUTINE STPNT

!----------------------------------------------------------------------
!----------------------------------------------------------------------

SUBROUTINE BCND(NDIM,PAR,ICP,NBC,U0,U1,FB,IJAC,DBC)
	!--------- ----

	! Boundary Conditions

	! Input arguments :
	!      NDIM   :   Dimension of the ODE system 
	!      PAR    :   Equation parameters
	!      ICP    :   Array indicating the free parameter(s)
	!      NBC    :   Number of boundary conditions
	!      U0     :   State variable values at the left boundary
	!      U1     :   State variable values at the right boundary
    
	! Values to be returned :
	!      FB     :   The values of the boundary condition functions 

	! Normally unused Jacobian arguments : IJAC, DBC (see manual)

	IMPLICIT NONE
	INTEGER, INTENT(IN) :: NDIM, ICP(*), NBC, IJAC
	DOUBLE PRECISION, INTENT(IN) :: PAR(*), U0(NDIM), U1(NDIM)
	DOUBLE PRECISION, INTENT(OUT) :: FB(NBC)
	DOUBLE PRECISION, INTENT(INOUT) :: DBC(NBC,*)

	!X FB(1)=
	!X FB(2)=

END SUBROUTINE BCND

!----------------------------------------------------------------------
!----------------------------------------------------------------------

SUBROUTINE ICND(NDIM,PAR,ICP,NINT,U,UOLD,UDOT,UPOLD,FI,IJAC,DINT)

	! Integral Conditions

	! Input arguments :
	!      NDIM   :   Dimension of the ODE system 
	!      PAR    :   Equation parameters
	!      ICP    :   Array indicating the free parameter(s)
	!      NINT   :   Number of integral conditions
	!      U      :   Value of the vector function U at `time' t

	! The following input arguments, which are normally not needed,
	! correspond to the preceding point on the solution branch
	!      UOLD   :   The state vector at 'time' t
	!      UDOT   :   Derivative of UOLD with respect to arclength
	!      UPOLD  :   Derivative of UOLD with respect to `time'

	! Normally unused Jacobian arguments : IJAC, DINT

	! Values to be returned :
	!      FI     :   The value of the vector integrand 

	IMPLICIT NONE
	INTEGER, INTENT(IN) :: NDIM, ICP(*), NINT, IJAC
	DOUBLE PRECISION, INTENT(IN) :: PAR(*)
	DOUBLE PRECISION, INTENT(IN) :: U(NDIM), UOLD(NDIM), UDOT(NDIM), UPOLD(NDIM)
	DOUBLE PRECISION, INTENT(OUT) :: FI(NINT)
	DOUBLE PRECISION, INTENT(INOUT) :: DINT(NINT,*)

END SUBROUTINE ICND

!----------------------------------------------------------------------
!----------------------------------------------------------------------


SUBROUTINE FOPT(NDIM,U,ICP,PAR,IJAC,FS,DFDU,DFDP)
	!--------- ----
	!
	! Defines the objective function for algebraic optimization problems
	!
	! Supplied variables :
	!      NDIM   :   Dimension of the state equation
	!      U      :   The state vector
	!      ICP    :   Indices of the control parameters
	!      PAR    :   The vector of control parameters
	!
	! Values to be returned :
	!      FS      :   The value of the objective function
	!
	! Normally unused Jacobian argument : IJAC, DFDP

	IMPLICIT NONE
	INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
	DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
	DOUBLE PRECISION, INTENT(OUT) :: FS
	DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM),DFDP(*)

END SUBROUTINE FOPT

!----------------------------------------------------------------------
!----------------------------------------------------------------------

SUBROUTINE PVLS(NDIM,U,PAR)
	!--------- ----

	IMPLICIT NONE
	INTEGER, INTENT(IN) :: NDIM
	DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM)
	DOUBLE PRECISION, INTENT(INOUT) :: PAR(*)
	DOUBLE PRECISION :: GETP,pi,realfm,imagfm,imagfm1
	DOUBLE PRECISION :: lw,lw1
	LOGICAL, SAVE :: first = .TRUE.
	DOUBLE PRECISION :: T
	INTEGER :: i

	!IF (first) THEN
		!CALL STPNT(NDIM,U,PAR,T)
		!first = .FALSE.
	!ENDIF

	PAR(25)=0.
	pi = 4*ATAN(1d0)
	i=1
	lw=100.
	lw1=101.
	DO WHILE(i < NDIM)
		realfm = GETP('EIG',I*2-1,U)
		IF (ABS(realfm) < lw) THEN
			lw = ABS(realfm)
			lw1 = ABS(GETP('EIG',(I+1)*2-1,U))
			imagfm1 = ABS(GETP('EIG',(I+1)*2,U))
			imagfm = ABS(GETP('EIG',I*2,U))
		END IF
		i=i+1
	END DO
	IF ((lw==lw1).AND.(imagfm1==imagfm).AND.(imagfm/=0.D0)) THEN
	PAR(25) = 2*pi/imagfm
	ENDIF
	!---------------------------------------------------------------------- 
	! NOTE : 
	! Parameters set in this subroutine should be considered as ``solution 
	! measures'' and be used for output purposes only.
	! 
	! They should never be used as `true'' continuation parameters. 
	!
	! They may, however, be added as ``over-specified parameters'' in the 
	! parameter list associated with the AUTO-Constant NICP, in order to 
	! print their values on the screen and in the ``p.xxx file.
	!
	! They may also appear in the list associated with AUTO-Constant NUZR.
	!
	!---------------------------------------------------------------------- 
	! For algebraic problems the argument U is, as usual, the state vector.
	! For differential equations the argument U represents the approximate 
	! solution on the entire interval [0,1]. In this case its values must 
	! be accessed indirectly by calls to GETP, as illustrated below.
	!---------------------------------------------------------------------- 
	!
	! Set PAR(2) equal to the L2-norm of U(1)
	!X PAR(2)=GETP('NRM',1,U)
	!
	! Set PAR(3) equal to the minimum of U(2)
	!X PAR(3)=GETP('MIN',2,U)
	!
	! Set PAR(4) equal to the value of U(2) at the left boundary.
	!X PAR(4)=GETP('BV0',2,U)
	!
	! Set PAR(5) equal to the pseudo-arclength step size used.
	!X PAR(5)=GETP('STP',1,U)
	!
	!---------------------------------------------------------------------- 
	! The first argument of GETP may be one of the following:
	!        'NRM' (L2-norm),     'MAX' (maximum),
	!        'INT' (integral),    'BV0 (left boundary value),
	!        'MIN' (minimum),     'BV1' (right boundary value).
	!
	! Also available are
	!   'STP' (Pseudo-arclength step size used).
	!   'FLD' (`Fold function', which vanishes at folds).
	!   'BIF' (`Bifurcation function', which vanishes at singular points).
	!   'HBF' (`Hopf function'; which vanishes at Hopf points).
	!   'SPB' ( Function which vanishes at secondary periodic bifurcations).
	!---------------------------------------------------------------------- 
END SUBROUTINE PVLS
