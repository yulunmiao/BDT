#include "make_BDT_var.h"
#define LUMI 59.74

int make_BDT_var(){
	//WWW
	make_BDT("WWW_2017.root","WWW_BDT_2017.root");
	make_BDT("WWW_2016.root","WWW_BDT_2016.root");
	make_BDT("WWW_filtered_2016.root","WWW_BDT_filtered_2016");
	make_BDT("WWW_filtered_2017.root","WWW_BDT_filtered_2017");
	//WWW
	return 0;
}
