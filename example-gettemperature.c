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


/// \file example-gettemperature.cpp
/// \brief This file contains code example for collecting temperature sensor data for Intel® Xeon Phi™ Coprocessor components.
/// \details This API is used to fetch different temperature sensor data and is available in both Intel® Xeon Phi™ Software Development Vehicle (formerly known as “Knights Ferry”) & KNC.
///  Intel® Xeon Phi™ Software Development Vehicle (formerly known as “Knights Ferry”) has 3 temperature sensors and Intel® Xeon Phi™ product codenamed “Knights Corner” has 8 temperature sensors.

///  API Used: MicGetTemperature
///  Header File for the API: MicThermalAPI.h
///  Data Structure Used: MicDeviceOnSystem from MicAccessTypes.h
///

#define MAX_DEVICES		(32)

#include <stdio.h>

#include "Types.h"
#include "MicAccessTypes.h"
#include "MicBasicTypes.h"
#include "MicAccessErrorTypes.h"
#include "MicAccessApi.h"
#include "MicThermalTypes.h"
#include "MicThermalAPI.h"	//API Specific Header File
#include "stdlib.h"

int main()
{
    MicDeviceOnSystem adaptersList[MAX_DEVICES];
    HANDLE accessHandle = NULL;
    U32 nAdapters = MAX_DEVICES;
    U32 adapterNum;
    U32 retVal = MIC_ACCESS_API_ERROR_UNKNOWN;
    U32 returningSize = 0;
    U32 *temperatures = (U32 *) NULL;
        int END_TEMP_SENSOR = eMicThermalVddq;
    int sensor;

    // Initialize the API class. Use eTARGET_SCIF_DRIVER target for Linux.
    retVal = MicInitAPI(&accessHandle, eTARGET_SCIF_DRIVER, adaptersList,
                &nAdapters);
    if (retVal != MIC_ACCESS_API_SUCCESS) {
        fprintf(stderr, "%s\n", MicGetErrorString(retVal));
        return retVal;
    }

    // To ensure that the adpterList doesn't overflow
    if (nAdapters < 0 || nAdapters >= MAX_DEVICES) {
        fprintf(stderr, "%s\n",
            MicGetErrorString(MIC_ACCESS_API_ERROR_UNKNOWN));
        return retVal;
    }
    // Iterate through the list of available cards (referred to as adapters
    // below), initialize them. Following this, call an API and then
    // close the adapter.
    for (adapterNum = 0; adapterNum < nAdapters; adapterNum++) {

        // Initialize adapter
        retVal =
            MicInitAdapter(&accessHandle, &adaptersList[adapterNum]);
        if (retVal != MIC_ACCESS_API_SUCCESS) {
            fprintf(stderr, "%s\n", MicGetErrorString(retVal));
            MicCloseAPI(&accessHandle);
            return retVal;
        }

                for (sensor = eMicThermalBoard; sensor <= END_TEMP_SENSOR;
             sensor++) {
            temperatures = (U32 *) malloc(returningSize);
            if (!temperatures)
            {
                MicCloseAdapter(&accessHandle);
                MicCloseAPI(&accessHandle);
                fprintf(stderr, "Error allocating temperature buffer\n");
                return MIC_ACCESS_API_INSUFFICIENT_MEMORY;
            }
            returningSize = 16;

            //API Invocation
            retVal =
                MicGetTemperature(accessHandle,
                                              (E_MIC_THERMAL_TYPE) sensor,
                          temperatures, &returningSize);
            if (retVal != MIC_ACCESS_API_SUCCESS) {
                fprintf(stderr, "%s\n",
                    MicGetErrorString(retVal));
                MicCloseAdapter(accessHandle);
                MicCloseAPI(&accessHandle);
                return retVal;
            }

            switch (sensor) {
                        case eMicThermalBoard:
                printf("Board Temperature: ");
                break;
                        case eMicThermalDevMem:
                                printf("DevMem Temperature: ");
                break;
                        case eMicThermalDie:
                printf("DIE Temperature: ");
                break;
                        case eMicThermalFin:
                printf("Inlet Temperature: ");
                break;
                        case eMicThermalFout:
                printf("Outlet Temperature: ");
                break;
                        case eMicThermalVccp:
                printf("Vccp Temperature: ");
                break;
                        case eMicThermalVddg:
                printf("Vddg Temperature: ");
                break;
                        case eMicThermalVddq:
                printf("Vddq Temperature: ");
                break;
            }

            printf("%u C\n", temperatures[0]);

            free(temperatures);
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
        fprintf(stderr, "%s\n", MicGetErrorString(retVal));
    }

    return retVal;
}

