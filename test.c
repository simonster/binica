#include <stdio.h>

main()
{
	double w[1000];
 	int i;

	for (i=0; i<1000; i++)
		w[i] = i;
       	
	printf("data length = %d\n",sizeof(w));
}

