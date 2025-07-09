
 Author        : Nino Jacobs, Catharina Hospital, Eindhoven, The Netherlands
 Created On    : 19-06-2025
 Software      : RayStation version 2024B

 Description   : GUI-based application for evaluating the impact of rectal air during
                 prostate cancer radiotherapy. The tool predicts whether the presence of
                 rectal air may lead to exceeding clinical dose thresholds.
                 If necessary, it advises intervention (e.g., patient removal and gas release)
                 before continuing treatment.

 Intended Use  : For patients treated under the in-house protocol PR66.
                 To be used when visible rectal air is observed on CBCT,
                 after a minimum of 5 delivered treatment fractions.

Files          : The Word document titled "" provides a description of how radiotherapy laborants can use the workflow.
                 It outlines the steps for initializing the workflow in RayStation and explains how the Python script (filename: "") is 
                 integrated into RayStation to determine whether intervention is necessary when rectal air is present.
