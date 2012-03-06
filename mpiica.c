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

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "mpiica.h"


/********************* Concatenate binary float matricies *********************/
/* Convert size doublereal elements of mat to floating point values and       */
/* appends values to the binary file specified by file.                       */
/*                                                                            */
/* file: FILE pointer (input)                                                 */
/* size: int (input)                                                          */
/* mat:  doublereal array [size] (input)                                      */

void fbc_matwrite(FILE *file, int size, doublereal *mat) {
	float *buffer;
	int i, items;

	if (!file) error("open failed");

	buffer = (float*)malloc(size*sizeof(float));

	for (i=0 ; i<size ; i++) buffer[i] = (float)mat[i];
	
	items = (int)fwrite(buffer,sizeof(float),size,file);
	if (items != size) error("invalid number of elements");
	
	free(buffer);
}


/********************* Concatenate ascii integer matricies ********************/
/* Append size integer values of mat as ascii numbers to the file specified   */
/* by file.                                                                   */
/*                                                                            */
/* file: FILE pointer (input)                                                 */
/* size: int (input)                                                          */
/* mat:  integer array [size] (input)                                         */

void iac_matwrite(FILE *file, int size, integer *mat) {
	int i;

	if (!file) error("open failed");
	for (i=0 ; i<size-1 ; i++) fprintf(file,"%d ",(int)(mat[i]));
	fprintf(file,"%d\n",(int)(mat[i]));
}


/************************** Send assignment to slave **************************/
/* Send assignment structure assign to slave process tid.                     */
/*                                                                            */
/* tid:    int (input)                                                        */
/* assign: tassign pointer (input)                                            */

void send_assign(int tid, tassign *assign) {
	int     bufid, datasize, i, toint[ASSIGN_NINT];
	double *tdbl, todbl[ASSIGN_NDBL];
	float  *tflt;

	toint[ASSIGN_ID]         = (int)(assign->id);
	toint[ASSIGN_CHANS]      = (int)(assign->chans);
	toint[ASSIGN_FRAMES]     = (int)(assign->frames);
	toint[ASSIGN_EPOCHS]     = (int)(assign->epochs);
	toint[ASSIGN_BIAS]       = (int)(assign->bias);
	toint[ASSIGN_SIGNS]      = (int)(assign->signs);
	toint[ASSIGN_EXTENDED]   = (int)(assign->extended);
	toint[ASSIGN_EXTBLOCKS]  = (int)(assign->extblocks);
	toint[ASSIGN_PDFSIZE]    = (int)(assign->pdfsize);
	toint[ASSIGN_NSUB]       = (int)(assign->nsub);
	toint[ASSIGN_VERBOSE]    = (int)(assign->verbose);
	toint[ASSIGN_BLOCK]      = (int)(assign->block);
	toint[ASSIGN_MAXSTEPS]   = (int)(assign->maxsteps);
	
	todbl[ASSIGN_LRATE]      = (double)(assign->lrate);
	todbl[ASSIGN_ANNEALSTEP] = (double)(assign->annealstep);
	todbl[ASSIGN_ANNEALDEG]  = (double)(assign->annealdeg);
	todbl[ASSIGN_NOCHANGE]   = (double)(assign->nochange);
	todbl[ASSIGN_MOMENTUM]   = (double)(assign->momentum);

	MPI_Send(toint,ASSIGN_NINT,MPI_INT,tid,1,MPI_COMM_WORLD);
	MPI_Send(todbl,ASSIGN_NDBL,MPI_DOUBLE,tid,1,MPI_COMM_WORLD);

	datasize = assign->chans * assign->frames * assign->epochs;
	if (sizeof(float) != sizeof(doublereal)) {
		tflt = (float*)malloc(datasize*sizeof(float));
		for (i=0 ; i<datasize ; i++) tflt[i] = (float)(assign->data[i]);
		MPI_Send(tflt,datasize,MPI_FLOAT,tid,1,MPI_COMM_WORLD);
		free(tflt);
	}
	else
		MPI_Send((float*)(assign->data),datasize,MPI_FLOAT,tid,1,MPI_COMM_WORLD);
		
	datasize = assign->chans * assign->chans;
	if (sizeof(double) != sizeof(doublereal)) {
		tdbl = (double*)malloc(datasize*sizeof(double));
		for (i=0 ; i<datasize ; i++) tdbl[i] = (double)(assign->weights[i]);
		MPI_Send(tdbl,datasize,MPI_DOUBLE,tid,1,MPI_COMM_WORLD);
		free(tdbl);
	}
	else
		MPI_Send((double*)(assign->weights),datasize,MPI_DOUBLE,tid,1,MPI_COMM_WORLD);
}


/*********************** Receive assignment from master ***********************/
/* Receive assignment structure assign from master.                           */
/*                                                                            */
/* assign: tassign pointer (output)                                           */

int receive_assign(tassign *assign) {
	int         bufid, tag, dummy, datasize, i, toint[ASSIGN_NINT];
	double     *tdbl, todbl[ASSIGN_NDBL];
	float      *tflt;
	MPI_Status  status;
	
	MPI_Recv(toint,ASSIGN_NINT,MPI_INT,-1,-1,MPI_COMM_WORLD,status);
	MPI_Recv(todbl,ASSIGN_NDBL,MPI_DOUBLE,-1,-1,MPI_COMM_WORLD,status);

	assign->id         = (integer)toint[ASSIGN_ID];
	assign->chans      = (integer)toint[ASSIGN_CHANS];
	assign->frames     = (integer)toint[ASSIGN_FRAMES];
	assign->epochs     = (integer)toint[ASSIGN_EPOCHS];
	assign->bias       = (integer)toint[ASSIGN_BIAS];
	assign->signs      = (integer)toint[ASSIGN_SIGNS];
	assign->extended   = (integer)toint[ASSIGN_EXTENDED];
	assign->extblocks  = (integer)toint[ASSIGN_EXTBLOCKS];
	assign->pdfsize    = (integer)toint[ASSIGN_PDFSIZE];
	assign->nsub       = (integer)toint[ASSIGN_NSUB];
	assign->verbose    = (integer)toint[ASSIGN_VERBOSE];
	assign->block      = (integer)toint[ASSIGN_BLOCK];
	assign->maxsteps   = (integer)toint[ASSIGN_MAXSTEPS];
	
	assign->lrate      = (doublereal)todbl[ASSIGN_LRATE];
	assign->annealstep = (doublereal)todbl[ASSIGN_ANNEALSTEP];
	assign->annealdeg  = (doublereal)todbl[ASSIGN_ANNEALDEG];
	assign->nochange   = (doublereal)todbl[ASSIGN_NOCHANGE];
	assign->momentum   = (doublereal)todbl[ASSIGN_MOMENTUM];

	datasize = assign->chans * assign->frames * assign->epochs;
	assign->data = (doublereal*)malloc(datasize*sizeof(doublereal));
	if (sizeof(float) != sizeof(doublereal)) {
		tflt = (float*)malloc(datasize*sizeof(float));
		MPI_Recv(tflt,datasize,MPI_FLOAT,-1,-1,MPI_COMM_WORLD,status);
		for (i=0 ; i<datasize ; i++) assign->data[i] = (doublereal)tflt[i];
		free(tflt);
	}
	else
		MPI_Recv((float*)(assign->data),datasize,MPI_FLOAT,-1,-1,MPI_COMM_WORLD,status);

	datasize = assign->chans * assign->chans;
	assign->weights = (doublereal*)malloc(datasize*sizeof(doublereal));
	if (sizeof(double) != sizeof(doublereal)) {
		tdbl = (double*)malloc(datasize*sizeof(double));
		MPI_Recv(tdbl,datasize,MPI_DOUBLE,-1,-1,MPI_COMM_WORLD,status);
		for (i=0 ; i<datasize ; i++) assign->weights[i] = (doublereal)tdbl[i];
		free(tdbl);
	}
	else
		MPI_Recv((double*)(assign->weights),datasize,MPI_DOUBLE,-1,-1,MPI_COMM_WORLD,status);

	return 1;
}


/*************************** Send results to master ***************************/
/* Send result structure result to master process tid.                        */
/*                                                                            */
/* tid:    int (input)                                                        */
/* result: tresult pointer (input)                                            */

void send_result(int tid, tresult *result) {
	int    *tint, bufid, datasize, i, toint[RESULT_NINT];
	double *tdbl, todbl[RESULT_NDBL];

	toint[RESULT_ID]    = (int)(result->id);
	toint[RESULT_CHANS] = (int)(result->chans);
	toint[RESULT_BIAS]  = (int)(result->bias != NULL);
	toint[RESULT_SIGNS] = (int)(result->signs != NULL);

	todbl[RESULT_LRATE] = (double)(result->lrate);
	
	MPI_Send(toint,RESULT_NINT,MPI_INT,tid,1,MPI_COMM_WORLD);
	MPI_Send(todbl,RESULT_NDBL,MPI_DOUBLE,tid,1,MPI_COMM_WORLD);

	datasize = result->chans * result->chans;
	if (sizeof(double) != sizeof(doublereal)) {
		tdbl = (double*)malloc(datasize*sizeof(double));
		for (i=0 ; i<datasize ; i++) tdbl[i] = (double)(result->weights[i]);
		MPI_Send(tdbl,datasize,MPI_DOUBLE,tid,1,MPI_COMM_WORLD);
		free(tdbl);
	}
	else
		MPI_Send((double*)(result->weights),datasize,MPI_DOUBLE,tid,1,MPI_COMM_WORLD);

	if (result->bias != NULL) {
		datasize = result->chans;
		if (sizeof(double) != sizeof(doublereal)) {
			tdbl = (double*)malloc(datasize*sizeof(double));
			for (i=0 ; i<datasize ; i++) tdbl[i] = (double)(result->bias[i]);
			MPI_Send(tdbl,datasize,MPI_DOUBLE,tid,1,MPI_COMM_WORLD);
			free(tdbl);
		}
		else
			MPI_Send((double*)(result->bias),datasize,MPI_DOUBLE,tid,1,MPI_COMM_WORLD);
	}

	if (result->signs != NULL) {
		datasize = result->chans;
		if (sizeof(int) != sizeof(integer)) {
			tint = (int*)malloc(datasize*sizeof(int));
			for (i=0 ; i<datasize ; i++) tint[i] = (int)(result->signs[i]);
			MPI_Send(tint,datasize,MPI_INT,tid,1,MPI_COMM_WORLD);
			free(tint);
		}
		else
			MPI_Send((int*)(result->signs),datasize,MPI_INT,tid,1,MPI_COMM_WORLD);
	}
}


/************************** Receive result from slave *************************/
/* Receive result structure result from the first slave to responding.        */
/*                                                                            */
/* result: tresult pointer (output)                                           */

int receive_result(tresult *result) {
	int    *tint, bufid, datasize, bias, signs, toint[RESULT_NINT];
	int     i, buflen, tag, tid;
	double *tdbl, todbl[RESULT_NDBL];

	MPI_Recv(toint,RESULT_NINT,MPI_INT,-1,-1,MPI_COMM_WORLD,status);
	MPI_Recv(todbl,RESULT_NDBL,MPI_DOUBLE,-1,-1,MPI_COMM_WORLD,status);

	result->id         = (integer)toint[RESULT_ID];
	result->chans      = (integer)toint[RESULT_CHANS];
	bias               = (integer)toint[RESULT_BIAS];
	signs              = (integer)toint[RESULT_SIGNS];

	result->lrate      = (doublereal)todbl[RESULT_LRATE];

	datasize = result->chans * result->chans;
	result->weights = (doublereal*)malloc(datasize*sizeof(doublereal));
	if (sizeof(double) != sizeof(doublereal)) {
		tdbl = (double*)malloc(datasize*sizeof(double));
		MPI_Recv(tdbl,datasize,MPI_DOUBLE,-1,-1,MPI_COMM_WORLD,status);
		for (i=0 ; i<datasize ; i++) result->weights[i] = (doublereal)tdbl[i];
		free(tdbl);
	}
	else
		MPI_Recv((double*)(result->weights),datasize,MPI_DOUBLE,-1,-1,MPI_COMM_WORLD,status);

	if (bias) {
		datasize = result->chans;
		result->bias = (doublereal*)malloc(datasize*sizeof(doublereal));
		if (sizeof(double) != sizeof(doublereal)) {
			tdbl = (double*)malloc(datasize*sizeof(double));
			MPI_Recv(tdbl,datasize,MPI_DOUBLE,-1,-1,MPI_COMM_WORLD,status);
			for (i=0 ; i<datasize ; i++) result->bias[i] = (doublereal)tdbl[i];
			free(tdbl);
		}
		else
			MPI_Recv((double*)(result->bias),datasize,MPI_DOUBLE,-1,-1,MPI_COMM_WORLD,status);
	}
	else
		result->bias = NULL;

	if (signs) {
		datasize = result->chans;
		result->signs = (integer*)malloc(datasize*sizeof(integer));
		if (sizeof(int) != sizeof(integer)) {
			tint = (int*)malloc(datasize*sizeof(int));
			MPI_Recv(tint,datasize,MPI_INT,-1,-1,MPI_COMM_WORLD,status);
			for (i=0 ; i<datasize ; i++) result->signs[i] = (integer)tint[i];
			free(tint);
		}
		else
			MPI_Recv((int*)(result->signs),datasize,MPI_INT,-1,-1,MPI_COMM_WORLD,status);
	}
	else
		result->signs = NULL;
	
	return tid;
}


/***************************** Main slave routine *****************************/
/* Repeatedly performs ICA on data segments. Assignments are received by the  */
/* process which decomposes transmitted data and returns the results to the   */
/* assigning master process. The sequence is repeated until a kill signal is  */
/* received.                                                                  */

int slave() {
	int         i, datasize, id;
	tassign     assign;
	tresult     result;
	doublereal *bias, *data, *weights;
	integer    *signs, chans, frames, epochs;

#ifndef CRAY
	setpriority(PRIO_PROCESS,0,SLAVE_NICE);
#endif

	while (receive_assign(&assign)) {
		id         = assign.id;
		data       = assign.data;
		weights    = assign.weights;
		chans      = assign.chans;
		frames     = assign.frames;
		epochs     = assign.epochs;

		extended   = assign.extended;
		extblocks  = assign.extblocks;
		pdfsize    = assign.pdfsize;
		nsub       = assign.nsub;
		verbose    = assign.verbose;
		block      = assign.block;
		maxsteps   = assign.maxsteps;
		lrate      = assign.lrate;
		annealstep = assign.annealstep;
		annealdeg  = assign.annealdeg;
		nochange   = assign.nochange;
		momentum   = assign.momentum;

		if (assign.bias)
			bias = (doublereal*)malloc(chans*sizeof(doublereal));
		else
			bias = NULL;

		if (assign.signs)
			signs = (integer*)malloc(chans*sizeof(integer));
		else
			signs = NULL;

		runica(data,weights,chans,frames,epochs,bias,signs);
	
		result.id      = id;
		result.weights = weights;
		result.chans   = chans;
		result.bias    = bias;
		result.signs   = signs;
		result.lrate   = lrate;

		send_result(0,&result);

		if (data != NULL) free(data);
		if (weights != NULL) free(weights);
		if (bias != NULL) free(bias);
		if (signs != NULL) free(signs);
	}
	return 0;
}


/*************************** Extracts a data segment **************************/
/* Exctracts the data segment corresponding to step from matrix data. The     */
/* resulting segment is stored in windata. Chans, frames and epochs denote    */
/* the dimensions of data. The window array specifies dimensions for grid     */
/* from which to extract windata                                              */
/*                                                                            */
/* data:     doublereal array [chans,frames*epoch] (input)                    */
/* datawin:  doublereal array -                                               */
/*                   [chans,window(FRAMEWINDOW)*window(EPOCHWINDOW)] (output) */
/* chans:    int (input)                                                      */
/* frames:   int (input)                                                      */
/* epochs:   int (input)                                                      */
/* step:     int (input)                                                      */
/* window:   int array [NWINDOW] (input)                                      */

doublereal *extract(doublereal *data, doublereal *windata, int chans, int frames, int epochs, int step, int *window) {
	int i, j, nep, ofs, size;

	nep = (epochs-window[EPOCHWINDOW])/window[EPOCHSTEP] + 1;
	ofs = frames * window[EPOCHSTEP] * (step%nep);
	ofs += window[FRAMESTEP] * (step/nep);
	size = window[FRAMEWINDOW]*chans;
	
	for (i=0,j=0 ; i<window[EPOCHWINDOW] ; i++,j+=size,ofs+=frames)
		memcpy(&windata[j],&data[ofs*chans],size*sizeof(doublereal));

	return windata;
}


/*************************** Extracts data baseline ***************************/
/* Exctracts and concatenates all baseline data contained in matrix data. The */
/* resulting data block is stored in basedata. Chans, frames and epochs       */
/* denote the dimensions of data. The window array contains, amongst others,  */
/* the baseline length.                                                       */
/*                                                                            */
/* data:     doublereal array [chans,frames*epoch] (input)                    */
/* basedata: doublereal array [chans,window(BASELINE)*epoch] (output)         */
/* chans:    int (input)                                                      */
/* frames:   int (input)                                                      */
/* epochs:   int (input)                                                      */
/* window:   int array [NWINDOW] (input)                                      */

doublereal *baseline(doublereal *data, doublereal *basedata, int chans, int frames, int epochs, int *window) {
	int i, j, k, size;

	size = window[BASELINE]*chans;
	for (i=0,j=0,k=0 ; i<epochs ; i++,j+=size,k+=frames*chans)
		memcpy(&basedata[j],&data[k],size*sizeof(doublereal));

	return basedata;
}


/************************** Executes ICA in parallel **************************/
/* A serial ICA decomposition is initially performed on the baseline portions */
/* of data, all succeeding decompositions are use this result as initial      */
/* weight estimate. A number of tasks are spawned on each node according to   */
/* host info speed values. Each node is repeatedly assigned jobs until all    */
/* decompositions in the frame/epoch grid specified by window are distributed.*/
/* Weights returned by slave processes are continously sorted by projected    */
/* variance stored and to the files specified by fnames. The following        */
/* externally accessible variables are assumed initialized: extended,         */
/* extblocks, pdfsize, nsub, verbose, block, maxsteps, lrate, annealstep,     */
/* annealdeg, nochange, and momentum. If the boolean variable extended is set,*/
/* signs must be defined (i.e. not NULL)                                      */
/*                                                                            */
/* data:    double array [ncomps,frames*epoch] (input)                        */
/* weights: double array [ncomps,ncomps] (input)                              */
/* sphere:  double array [ncomps,ncomps] (input)                              */
/* eigv:    double array [ncomps,chans] (input)                               */
/* chans:   integer (input)                                                   */
/* ncomps:  integer (input)                                                   */
/* frames:  integer (input)                                                   */
/* epochs:  integer (input)                                                   */
/* window:  integer array [NWINDOW] (input)                                   */
/* bias:    double array [ncomps] (dummy) or NULL                             */
/* signs:   integer array [ncomps] (dummy) or NULL                            */
/* fnames:  char array [3] of array                                           */

void mpiica(doublereal *data, doublereal *weights, doublereal *sphere, doublereal *eigv, integer chans, integer ncomps, integer frames, integer epochs, int *window, doublereal *bias, integer *signs, char **fnames) {
	int                 i, j, datasize, maxep, maxfr, segs, speed, spwnd = 0, id = 0;
	int                 tid, nproc, ntask, last = 0;
	char               *name;
	FILE               *fids[3];
	integer           **srec;
	doublereal        **wrec, **brec, *basedata, *windata, *prjdata, deflr = lrate;
	tresult             result;
	tassign             assign;

	MPI_Comm_size(MPI_COMM_WORLD,&nproc);

	basedata = (doublereal*)malloc(ncomps*epochs*window[BASELINE]*sizeof(doublereal));
	baseline(data,basedata,ncomps,frames,epochs,window);
	runica(basedata,weights,ncomps,1,window[BASELINE]*epochs,bias,signs);
	free(basedata);

	assign.weights = weights;
	assign.chans = ncomps;
	assign.frames = window[FRAMEWINDOW];
	assign.epochs = window[EPOCHWINDOW];
	assign.bias = (int)(bias!=NULL);
	assign.signs = (int)(signs!=NULL);
	assign.extended = extended;
	assign.extblocks = extblocks;
	assign.pdfsize = pdfsize;
	assign.nsub = nsub;
	assign.verbose = verbose;
	assign.block = block;
	assign.maxsteps = maxsteps;
	assign.lrate = deflr;
	assign.annealstep = annealstep;
	assign.annealdeg = annealdeg;
	assign.nochange = nochange;
	assign.momentum = momentum;
	
	datasize = window[FRAMEWINDOW] * window[EPOCHWINDOW] * ncomps;
	windata = (doublereal*)malloc(datasize*sizeof(doublereal));
	prjdata = (doublereal*)malloc(datasize*sizeof(doublereal));
	
	maxep = (epochs-window[EPOCHWINDOW])/window[EPOCHSTEP] + 1;
	maxfr = (frames-window[FRAMEWINDOW])/window[FRAMESTEP] + 1;
	segs = maxep*maxfr;

	for (i=1 ; i<nproc ; i++,id++) {
		assign.id   = id;
		assign.data = extract(data,windata,(int)ncomps,(int)frames,(int)epochs,(int)id,window);
		send_assign(i,&assign);
	}
	ntask = nproc-1;

	wrec = (doublereal**)malloc(segs*sizeof(doublereal*));
	brec = (doublereal**)malloc(segs*sizeof(doublereal*));
	srec = (integer**)malloc(segs*sizeof(integer*));
	for (i=0 ; i<segs ; i++) {
		wrec[i] = NULL;
		brec[i] = NULL;
		srec[i] = NULL;
	}

	for (i=0 ; i<3 ; i++) fids[i] = NULL;
	if (fnames[0] != NULL) fids[0] = fopen(fnames[0],"wb");
	if (fnames[1] != NULL) fids[1] = fopen(fnames[1],"wb");
	if (fnames[2] != NULL) fids[2] = fopen(fnames[2],"wt");

	for (i=0 ; i<segs ; i++) {
		tid = receive_result(&result);
		
		printf("Received id %d\n",result.id);
		wrec[result.id] = result.weights;
		brec[result.id] = result.bias;
		srec[result.id] = result.signs;

		if (ntask < segs) {
			assign.id = id;
			assign.data = extract(data,windata,(int)ncomps,(int)frames,(int)epochs,(int)id,window);
			send_assign(tid,&assign);
			ntask++;
			id++;
		}

		datasize = window[FRAMEWINDOW] * window[EPOCHWINDOW];
		extract(data,windata,(int)ncomps,(int)frames,(int)epochs,(int)result.id,window);
		geproj(windata,result.weights,(integer)ncomps,(integer)datasize,prjdata);
		if (eigv)
			varsort(prjdata,result.weights,sphere,&eigv[chans*(chans-ncomps)],result.bias,result.signs,(integer)ncomps,(integer)datasize,(integer)chans);
		else
			varsort(prjdata,result.weights,sphere,NULL,result.bias,result.signs,(integer)ncomps,(integer)datasize,(integer)chans);

		while (last<segs && wrec[last]!=NULL) {
			if (fids[0]!=NULL && wrec[last]!=NULL) fbc_matwrite(fids[0],chans*ncomps,wrec[last]);
			if (fids[1]!=NULL && brec[last]!=NULL) fbc_matwrite(fids[1],ncomps,brec[last]);
			if (fids[2]!=NULL && srec[last]!=NULL) iac_matwrite(fids[2],ncomps,srec[last]);
			if (wrec[last] != NULL) free(wrec[last]);
			if (brec[last] != NULL) free(brec[last]);
			if (srec[last] != NULL) free(srec[last]);
			last++;
		}
	}

	for (i=0 ; i<3 ; i++) 
		if (fids[i] != NULL) fclose(fids[i]);

	if (wrec != NULL) free(wrec);
	if (brec != NULL) free(brec);
	if (srec != NULL) free(srec);
	if (windata != NULL) free(windata);
	if (prjdata != NULL) free(prjdata);
}
