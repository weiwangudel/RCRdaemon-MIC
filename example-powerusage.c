/*
* Copyright (C) Intel Corporation (2012)
*
* This file is subject to the Intel Sample Source Code License. 
* A copy of the Intel Sample Source Code License is included. 
*
* Linux OS: 	/opt/intel/mic/LICENSE/ 
* Windows: 	C:\Program Files\Intel\MPSS\
*
*/


/// \file example-powerusage.cpp
/// \brief This file contains code example for retrieving power readings of Intel® Xeon Phi™ Coprocessors and components.
/// \details API Used: MicGetPowerUsage
///  Header File for the API: MicPowerManagerAPI.h
///  Data Structures Used: MicDeviceOnSystem from MicAccessTypes.h
///                        MicPwrUsage from MicPowerManager.h

#define MAX_DEVICES		(32)

#include <stdio.h>

#include "Types.h"
#include "MicAccessTypes.h"
#include "MicBasicTypes.h"
#include "MicAccessErrorTypes.h"
#include "MicAccessApi.h"

#include "MicPowerManagerAPI.h"

int main()
{
	MicDeviceOnSystem adaptersList[MAX_DEVICES];
	HANDLE accessHandle = NULL;
	U32 nAdapters = MAX_DEVICES;
	U32 adapterNum = 0;
	U32 retVal = MIC_ACCESS_API_ERROR_UNKNOWN;
	MicPwrUsage powerUsage;

	// Initialize the API class. Use eTARGET_SCIF_DRIVER target for Linux.
	retVal = MicInitAPI(&accessHandle, eTARGET_SCIF_DRIVER, adaptersList,
			    &nAdapters);
	if (retVal != MIC_ACCESS_API_SUCCESS) {
		fprintf(stderr, "%s\n", MicGetErrorString(retVal));
		MicCloseAPI(&accessHandle);
		return retVal;
	}
	//Ensure that adapterList doesn't overflow
	if (nAdapters < 0 || nAdapters >= MAX_DEVICES) {
		fprintf(stderr, "%s\n",
			MicGetErrorString(MIC_ACCESS_API_ERROR_UNKNOWN));
		MicCloseAPI(&accessHandle);
		return retVal;
	}
	// Iterate through the list of available cards (referred to as adapters
	// below), initialize them. Following this, call the API to get uOS
	// version and then close the adapter.
	for (adapterNum = 0; adapterNum < nAdapters; adapterNum++) {
		// Initialize adapter
		retVal =
		    MicInitAdapter(&accessHandle, &adaptersList[adapterNum]);
		if (retVal != MIC_ACCESS_API_SUCCESS) {
			MicCloseAPI(&accessHandle);
			fprintf(stderr, "%s\n", MicGetErrorString(retVal));
			return retVal;
		}
        while (1) {
		// API call example: get and display the power usage.
		retVal = MicGetPowerUsage(accessHandle, &powerUsage);
		if (retVal != MIC_ACCESS_API_SUCCESS) {
			fprintf(stderr, "%s\n", MicGetErrorString(retVal));
			MicCloseAdapter(accessHandle);
			MicCloseAPI(&accessHandle);
			return retVal;
		}
		printf("Current Power Usage: %u\n", powerUsage.total0.prr);
        }
		// Close adapter
		retVal = MicCloseAdapter(accessHandle);
		if (retVal != MIC_ACCESS_API_SUCCESS) {
			fprintf(stderr, "%s\n", MicGetErrorString(retVal));
			MicCloseAPI(&accessHandle);
			return retVal;
		}
	}

	retVal = MicCloseAPI(&accessHandle);
	if (retVal != MIC_ACCESS_API_SUCCESS) {
	}
	return retVal;
}

