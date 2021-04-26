#include <cstdio>
#include "pwr.h" // link with -lpwr

int main()
{
	PWR_Cntxt cntxt = nullptr;
	PWR_Obj obj = nullptr;
	double clk0 = 0.0;
	PWR_Time ts0 = 0;

	PWR_CntxtInit(PWR_CNTXT_DEFAULT, PWR_ROLE_APP, "app", &cntxt);

	PWR_CntxtGetObjByName(cntxt, "plat.node.cpu", &obj);

	PWR_ObjAttrGetValue(obj, PWR_ATTR_FREQ, &clk0, &ts0);

	printf("%f MHz\n", clk0*1.e-6);

	PWR_CntxtDestroy(cntxt);
	return 0;
}
