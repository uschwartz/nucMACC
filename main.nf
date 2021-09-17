#!/usr/bin/env nextflow

 /*
 ===============================================================================
                      nextflow based nucMACC pipeline
 ===============================================================================
Authors:
Uwe Schwartz <uwe.schwartz@ur.de>
 -------------------------------------------------------------------------------
 */

nextflow.enable.dsl = 2

 //                           show settings

 if (!params.help) {
         include{settings} from './modules/setting'
         settings()
 }


 //                       help message



 // Show help message
 if (params.help) {
     include{helpMessage} from './modules/help'
     helpMessage()
     exit 0
 }

 //                      workflow
