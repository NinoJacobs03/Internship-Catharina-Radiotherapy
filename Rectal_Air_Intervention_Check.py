# *************************************************************************************************
# Script Name   : Workflow_rectal_air_CBCT
# Author        : Nino Jacobs, Catharina Hospital, Eindhoven, The Netherlands
# Created On    : 19-06-2025
# Software      : RayStation version 2024B
#
# Description   : GUI-based application for evaluating the impact of rectal air during
#                 prostate cancer radiotherapy. The tool predicts whether the presence of
#                 rectal air may lead to exceeding clinical dose thresholds.
#                 If necessary, it advises intervention (e.g., patient removal and gas release)
#                 before continuing treatment.
#
# Intended Use  : For patients treated under the in-house protocol PR66.
#                 To be used when visible rectal air is observed on CBCT,
#                 after a minimum of 5 delivered treatment fractions.
# *************************************************************************************************

from connect import *
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
from System.Windows import MessageBox
import tkinter as tk
from tkinter import ttk, font
import sys

# Append path to access general RayStation libraries
sys.path.append(sScriptsDir + r'\General\Libraries')
from RayScriptModule import showmessage


class Model:
    def __init__(self):
        # Retrieve current plan, beam set, exam, case and patient from RayStation
        self.plan = get_current('Plan')
        self.beam_set = get_current('BeamSet')
        self.exam = self.beam_set.GetPlanningExamination()
        self.case = get_current('Case')
        self.patient = get_current('Patient')

    def get_fraction_info(self):
        # Extract and format patient's full name
        name = self.patient.GetAlphabeticPatientName()
        full_name = f"{name['FirstName']} {name['LastName']}"

        # Count total and delivered treatment fractions
        total_fractions = len(self.plan.TreatmentCourse.TreatmentFractions)
        delivered_fractions = sum(
            1 for f in self.case.TreatmentDelivery.TreatmentCourse.TreatmentFractions if f.Status == 'Delivered'
        )

        # Subtract 1 to exclude the undelivered fraction with air (falsely marked as 'Delivered' in RayStation)
        return full_name, total_fractions, delivered_fractions - 1

    def evaluate_rectal_air(self, total_fractions, delivered_fractions):
        # Require minimum number of delivered fractions for reliable prediction
        if delivered_fractions < 5:
            MessageBox.Show(
                f"Only {delivered_fractions} fractions are delivered.\n\n"
                "A minimum of 5 delivered fractions is required for a reliable prediction.\n"
                "Intervention is necessary: the patient must be taken off the treatment table and\n"
                "an attempt should be made to remove rectal air before proceeding.\n\n"
                "If there are further uncertainties, consult the responsible physician.",
                "Insufficient Fractions – Intervention Required"
            )
            exit()

        roi_name = 'Rectum'
        threshold_v6000 = 6000  # Threshold dose for V6000 (cGy)
        threshold_v6200 = 6200  # Threshold dose for V6200 (cGy)
        volume_limit_v6200 = 2.0  # Clinical goal: max 2 cc at 6200 cGy
        dose_limit_v6000 = 3.0   # Clinical goal: max 3% volume at 6000 cGy

        CG_v6000 = []
        CG_v6200 = []

        # Retrieve dose grid and volume info for ROI (Rectum)
        dose_grid = self.beam_set.FractionDose.GetDoseGridRoi(RoiName=roi_name)
        voxel_indices = dose_grid.RoiVolumeDistribution.VoxelIndices
        relative_volumes = dose_grid.RoiVolumeDistribution.RelativeVolumes

        # Get shape of dose grid based on planning CT (deformed dose grid)
        grid_shape = self.case.TreatmentDelivery.TreatmentCourse.TreatmentFractions[0] \
            .EstimatedFractionDoseOnTotalDoseExamination.InDoseGrid.NrVoxels

        # Calculate voxel volume in cc
        voxel_volume_cc = (
            self.case.TreatmentDelivery.TreatmentCourse.TreatmentFractions[0]
            .EstimatedFractionDoseOnTotalDoseExamination.InDoseGrid.VoxelSize.x *
            self.case.TreatmentDelivery.TreatmentCourse.TreatmentFractions[0]
            .EstimatedFractionDoseOnTotalDoseExamination.InDoseGrid.VoxelSize.y *
            self.case.TreatmentDelivery.TreatmentCourse.TreatmentFractions[0]
            .EstimatedFractionDoseOnTotalDoseExamination.InDoseGrid.VoxelSize.z
        )

        # Initialize 3D matrix for cumulative delivered dose
        delivered_dose = np.zeros([grid_shape.z, grid_shape.y, grid_shape.x])

        # Accumulate dose from all delivered fractions (excluding undelivered 'air' fraction)
        for i in range(delivered_fractions):
            delivered_dose += self.case.TreatmentDelivery.TreatmentCourse.TreatmentFractions[i] \
                .EstimatedFractionDoseOnTotalDoseExamination.DoseValues.DoseData

        # Calculate mean delivered fractions dose and predicted total dose
        mean_dose = delivered_dose / delivered_fractions
        accumulated_dose = mean_dose * total_fractions

        # Extract doses within the ROI (rectum)
        roi_doses = accumulated_dose.flatten()[voxel_indices]

        # Calculate projected dose metrics
        current_v6000 = np.sum(relative_volumes[roi_doses >= threshold_v6000]) * 100   # V6000 in %
        current_v6200 = np.sum(roi_doses >= threshold_v6200) * voxel_volume_cc         # V6200 in cc

        # If any clinical goal is exceeded, notify and exit
        if current_v6000 >= dose_limit_v6000 or current_v6200 >= volume_limit_v6200:
            MessageBox.Show(
                f"Predicted Rectal V6000 cGy = {current_v6000:.2f}%\n"
                f"Predicted Rectal V6200 cGy = {current_v6200:.2f} cc\n\n"
                "One or both clinical goals are exceeded.\n"
                "Intervention is necessary: the patient must be taken off the treatment table and "
                "an attempt should be made to remove rectal air before proceeding.\n\n"
                "If there are further uncertainties, consult the responsible physician.",
                "Rectal Dose Too High – Intervention Required"
            )
            exit()

        # Get dose distribution for fraction with rectal air (currently not delivered)
        air_dose = self.case.TreatmentDelivery.TreatmentCourse.TreatmentFractions[delivered_fractions] \
            .EstimatedFractionDoseOnTotalDoseExamination.DoseValues.DoseData

        remaining_fractions = total_fractions - delivered_fractions

        # Simulate impact of increasing number of remaining fractions with rectal air
        for x in range(1, remaining_fractions + 1):
            y = total_fractions - x
            future_total_dose = y * mean_dose + x * air_dose
            simulated_doses = future_total_dose.flatten()[voxel_indices]

            v6000 = np.sum(relative_volumes[simulated_doses >= threshold_v6000]) * 100
            v6200 = np.sum(simulated_doses >= threshold_v6200) * voxel_volume_cc

            CG_v6000.append(v6000)
            CG_v6200.append(v6200)

        self.plot_prediction(CG_v6000, CG_v6200, delivered_fractions, total_fractions)

    def plot_prediction(self, CG_v6000, CG_v6200, delivered_fractions, total_fractions):

        # Inform user that results should be reviewed by the physician
        MessageBox.Show(
            "The following two graphs represent the predicted rectal dose levels:\n"
            "- V6000 cGy (%)\n"
            "- V6200 cGy (cc)\n\n"
            "These plots must be shown to the physician to determine appropriate clinical action.",
            "Physician Review Required"
        )

        fraction_numbers = list(range(delivered_fractions + 1, total_fractions + 1))
        fraction_numbers_air = list(range(1, total_fractions - delivered_fractions + 1))

        # Plot predicted V6000 cGy (%)
        fig1, ax1 = plt.subplots(figsize=(10, 6))
        ax1.set_xticks(fraction_numbers)
        ax1.plot(fraction_numbers, CG_v6000, marker='o', linestyle='-', color='black', label="V6000 cGy (%)")
        ax1.axhline(y=3, color='r', linestyle='--', label="Clinical goal threshold (3%)")
        ax1.set_xlabel("Fraction Number", fontsize=12)
        ax1.set_ylabel("Rectal V6000 cGy (%)", fontsize=12)
        ax1.grid(True)
        ax1.legend()

        # Add secondary x-axis showing number of air-containing fractions
        ax2 = ax1.twiny()
        ax2.set_xlim(ax1.get_xlim())
        ax2.set_xticks(fraction_numbers)
        ax2.set_xticklabels(fraction_numbers_air)
        ax2.set_xlabel("Number of rectal air fractions substituted", fontsize=12)

        plt.tight_layout()
        plt.show()

        # Plot predicted V6200 cGy (cc)
        fig2, ax3 = plt.subplots(figsize=(10, 6))
        ax3.set_xticks(fraction_numbers)
        ax3.plot(fraction_numbers, CG_v6200, marker='o', linestyle='-', color='black', label="V6200 cGy (cc)")
        ax3.axhline(y=2, color='r', linestyle='--', label="Clinical goal threshold (2 cc)")
        ax3.set_xlabel("Fraction Number", fontsize=12)
        ax3.set_ylabel("Rectal V6200 cGy (cc)", fontsize=12)
        ax3.grid(True)
        ax3.legend()

        # Add secondary x-axis showing number of air-containing fractions
        ax4 = ax3.twiny()
        ax4.set_xlim(ax3.get_xlim())
        ax4.set_xticks(fraction_numbers)
        ax4.set_xticklabels(fraction_numbers_air)
        ax4.set_xlabel("Number of rectal air fractions substituted", fontsize=12)

        plt.tight_layout()
        plt.show()


class View(ttk.Frame):
    def __init__(self, parent):
        super().__init__(parent)
        self.parent = parent
        self.bold_font = font.Font(family="Segoe UI", size=11, weight="bold")
        self._create_widgets()

    def _create_widgets(self):
        row = 0

        # Display plan name
        ttk.Label(self, text="Plan Name:", font=self.bold_font).grid(row=row, column=0, sticky=tk.W, padx=5, pady=3)
        ttk.Label(self, text=get_current('Plan').Name).grid(row=row, column=1, sticky=tk.W + tk.E, padx=5, pady=3)
        row += 1

        # Separator
        ttk.Separator(self, orient='horizontal').grid(row=row, column=0, columnspan=2, sticky=tk.EW, pady=(5, 10))
        row += 1

        # Patient and treatment info
        ttk.Label(self, text="Patient Name:", font=self.bold_font).grid(row=row, column=0, sticky=tk.W, padx=5, pady=3)
        self.patient_name_label = ttk.Label(self, text="")
        self.patient_name_label.grid(row=row, column=1, sticky=tk.W + tk.E, padx=5, pady=3)
        row += 1

        ttk.Label(self, text="Total Treatment Fractions:", font=self.bold_font).grid(row=row, column=0, sticky=tk.W, padx=5, pady=3)
        self.total_fractions_label = ttk.Label(self, text="")
        self.total_fractions_label.grid(row=row, column=1, sticky=tk.W + tk.E, padx=5, pady=3)
        row += 1

        ttk.Label(self, text="Fractions Delivered:", font=self.bold_font).grid(row=row, column=0, sticky=tk.W, padx=5, pady=3)
        self.delivered_fractions_label = ttk.Label(self, text="")
        self.delivered_fractions_label.grid(row=row, column=1, sticky=tk.W + tk.E, padx=5, pady=3)
        row += 1

        # Analysis header
        ttk.Separator(self, orient='horizontal').grid(row=row, column=0, columnspan=2, sticky=tk.EW, pady=(10, 10))
        row += 1

        ttk.Label(self, text="Analysis:", font=self.bold_font).grid(row=row, column=0, sticky=tk.W, padx=5, pady=3)
        ttk.Label(self, text="Evaluate rectal air impact").grid(row=row, column=1, sticky=tk.W + tk.E, padx=5, pady=3)
        row += 1

        # Goals evaluated
        ttk.Label(self, text="Clinical Goals Evaluated:", font=self.bold_font).grid(row=row, column=0, sticky=tk.W, padx=5, pady=3)
        ttk.Label(self, text="V6000 cGy (%) and V6200 cGy (cc)").grid(row=row, column=1, sticky=tk.W + tk.E, padx=5, pady=3)
        row += 1

        # Buttons
        ttk.Separator(self, orient='horizontal').grid(row=row, column=0, columnspan=2, sticky=tk.EW, pady=(10, 10))
        row += 1

        button_frame = ttk.Frame(self)
        button_frame.grid(row=row, column=0, columnspan=2, pady=(5, 5), sticky=tk.E)

        ttk.Button(button_frame, text="Cancel", command=self._cancel_clicked).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Analyse", command=self._analyse_clicked).pack(side=tk.LEFT, padx=5)

    def set_controller(self, controller):
        self.controller = controller

    def initialize_data(self):
        # Load and display patient data
        name, total, delivered = self.controller.get_fraction_info()
        self.patient_name_label.config(text=name)
        self.total_fractions_label.config(text=total)
        self.delivered_fractions_label.config(text=delivered)

    def _analyse_clicked(self):
        # Start analysis on click
        self.controller.run_analysis(
            int(self.total_fractions_label.cget("text")),
            int(self.delivered_fractions_label.cget("text"))
        )

    def _cancel_clicked(self):
        # Close application
        self.parent.destroy()


class Controller:
    def __init__(self, model, view):
        self.model = model
        self.view = view

    def get_fraction_info(self):
        return self.model.get_fraction_info()

    def run_analysis(self, total_fractions, delivered_fractions):
        self.model.evaluate_rectal_air(total_fractions, delivered_fractions)


class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title('Rectal Air Analysis')
        self.resizable(False, False)

        # Apply custom RayStation theme
        s = ttk.Style()
        self.tk.call('source', sScriptsDir + '/RayStationTheme/awdark.tcl')
        s.theme_use('awdark')

        model = Model()
        view = View(self)
        view.grid(row=0, column=0, padx=10, pady=10)

        controller = Controller(model, view)
        view.set_controller(controller)
        view.initialize_data()


if __name__ == "__main__":
    app = App()
    app.mainloop()
