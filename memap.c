/******************************************************************************/
/* For all support and information please contact:                            */
/*                                                                            */
/*   Sigurd Enghoff                                                           */
/*   The Salk Institute, CNL                                                  */
/*   enghoff@salk.edu                                                         */
/*                                                                            */
/* Additional ICA software:                                                   */
/*   http://www.cnl.salk.edu/~enghoff/                                        */
/*                                                                            */
/******************************************************************************/

#include <sys/types.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stddef.h>
#include <stdio.h>

#ifndef ALPHA

void *mapmalloc(int size) {
	int     fd, i;
	caddr_t base;
	char    dummy = 0;
	
	/*fd = open(tmpnam(NULL),O_RDWR|O_CREAT|O_EXCL);
	for (i=0 ; i<size ; i++) write(fd,&dummy,1);
	lseek(fd,SEEK_SET,0);*/
	
	fd = open("/dev/zero",O_RDWR);
	base = mmap(NULL,(size_t)(size),PROT_READ|PROT_WRITE,MAP_PRIVATE,fd,0);
	close(fd);
	return (void*)base;
}

#else

void *mapmalloc(int size) {
	caddr_t base;
	
	base = mmap(NULL,(size_t)(size),PROT_READ|PROT_WRITE,MAP_ANONYMOUS|MAP_VARIABLE,-1,0);
	return (void*)base;
}

#endif

void mapfree(void *addr, int size) {
	munmap((caddr_t)addr,size);
}

