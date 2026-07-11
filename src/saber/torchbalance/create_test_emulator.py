#!/usr/bin/env python3
"""
Create TorchScript emulators for testing TorchBalance.

This script demonstrates a simple linear emulator for testing the TorchBalance operator.
The model is intentionally simple (linear) for predictable, verifiable behavior in ctests.

For real ML models, Jacobians would be computed via:
- Automatic differentiation during training (PyTorch autograd)
- Pre-computed and stored in the model
- Computed on-demand using finite differences or analytical derivatives

This example shows the structure needed for TorchBalance compatibility.

Usage:
    python create_test_emulator.py  # Create test suite
"""

import torch
import torch.nn as nn
import argparse
from pathlib import Path
from typing import List


class SimpleLinearEmulator(nn.Module):
    """
    Simple linear emulator for testing TorchBalance.

    This model computes:
        output = a * input1 + b * input2 + c

    Where a, b, c are fixed constants. This demonstrates the structure needed
    for TorchBalance compatibility. Real ML emulators would have more complex
    forward passes and would compute Jacobians via automatic differentiation
    during training, storing the necessary information for efficient inference.
    """

    def __init__(self, num_inputs=2, jacobian_values=None):
        """
        Initialize the emulator.

        Args:
            num_inputs: Number of input variables
            jacobian_values: List of constant values for linear coefficients (one per input)
        """
        super().__init__()

        if jacobian_values is None:
            jacobian_values = [1.0] * num_inputs

        # Store linear coefficients as a single tensor (not trainable)
        self.register_buffer('coefficients', torch.tensor(jacobian_values))

        # Bias term (fixed)
        self.register_buffer('bias', torch.tensor(0.0))

        self.num_inputs = num_inputs

    def forward(self, inputs):
        """
        Forward pass: linear combination of inputs.

        Note: Required for TorchScript tracing but NOT called by C++ code.
        Only jac_physical() is called by C++.

        Args:
            inputs: Concatenated input tensor [batch_size, num_inputs]

        Returns:
            Output tensor [batch_size, 1]
        """
        # Linear combination: sum(coeff_i * input_i) + bias
        # inputs is [batch_size, num_inputs], coefficients is [num_inputs]
        output = torch.matmul(inputs, self.coefficients.unsqueeze(1)) + self.bias
        return output

    @torch.jit.export
    def jac_physical(self, inputs, mask):
        """
        Compute Jacobians with respect to all inputs in physical units.

        This is the ONLY method called by the C++ TorchBalanceEmulator code.

        For this simple linear model, the Jacobians are just the coefficients.
        In a real ML model, these would be computed via automatic differentiation
        using torch.autograd.grad() during model training/setup, then stored or
        computed on-the-fly.

        Note: Full autograd within TorchScript has limitations. For production ML
        emulators, Jacobians are typically:
        1. Pre-computed during training and stored
        2. Computed using finite differences
        3. Computed via custom analytical derivatives
        4. Or computed in Python and the results traced into TorchScript

        Args:
            inputs: Input tensor with shape [batch_size, num_inputs]
            mask: Mask tensor with shape [batch_size, 1]. Jacobian will be multiplied
                  by the mask (zeros where masked). Pass a tensor of ones for no masking.

        Returns:
            Jacobian tensor with shape [batch_size, 1, num_inputs]
            where jac[i, 0, j] = d(output)/d(input_j) for sample i
        """
        batch_size = inputs.shape[0]

        # For this linear model: output = coefficients · inputs + bias
        # Therefore: d(output)/d(input_j) = coefficient_j

        # Create Jacobian tensor: [batch_size, 1 output, num_inputs]
        jac = torch.zeros(batch_size, 1, self.num_inputs)

        # Fill with coefficients (constant for all samples in linear model)
        # Shape: coefficients is [num_inputs], we expand to [batch_size, num_inputs]
        jac[:, 0, :] = self.coefficients.unsqueeze(0).expand(batch_size, -1)

        # Apply mask (always provided)
        # Expand mask from [batch_size, 1] to [batch_size, 1, num_inputs]
        mask_expanded = mask.unsqueeze(2).expand(-1, 1, self.num_inputs)
        jac = jac * mask_expanded

        return jac


def create_emulator_with_metadata(
    output_var_name,
    input_var_names,
    jacobian_values=None,
    input_level=0,
    output_level=0,
    output_path="emulator.ts"
):
    """
    Create a TorchScript emulator with metadata for TorchBalance.

    Args:
        output_var_name: Name of output variable (e.g., "air_temperature")
        input_var_names: List of input variable names (e.g., ["sea_water_potential_temperature", "water_vapor_mixing_ratio_wrt_moist_air"])
        jacobian_values: List of constant Jacobian values (one per input)
        input_level: Level index for input variables (default: 0 for surface)
        output_level: Level index for output variable (default: 0 for surface)
        output_path: Path to save the TorchScript model
    """
    num_inputs = len(input_var_names)

    if jacobian_values is None:
        # Default: use different values for each Jacobian for easier testing
        jacobian_values = [0.5 + 0.1 * i for i in range(num_inputs)]

    # Create the model
    model = SimpleLinearEmulator(num_inputs=num_inputs, jacobian_values=jacobian_values)
    model.eval()

    # Create example inputs for tracing (batch_size=2)
    example_input = torch.randn(2, num_inputs)

    # Trace the model
    traced_model = torch.jit.trace(model, example_input)

    # Create a wrapper to hold both the model and metadata
    # This ensures metadata is properly serialized in TorchScript
    class EmulatorWithMetadata(torch.nn.Module):
        def __init__(self, traced_model, input_names, output_names, input_levels, output_levels):
            super().__init__()
            self.model = traced_model
            # Store metadata as constants (will be serialized)
            self.input_names = torch.jit.Attribute(input_names, List[str])
            self.output_names = torch.jit.Attribute(output_names, List[str])
            self.input_levels = torch.jit.Attribute(input_levels, List[int])
            self.output_levels = torch.jit.Attribute(output_levels, List[int])

        def forward(self, x):
            return self.model(x)

        @torch.jit.export
        def jac_physical(self, x, mask):
            return self.model.jac_physical(x, mask)

    # Create wrapper with metadata
    wrapper = EmulatorWithMetadata(
        traced_model,
        input_names=input_var_names,
        output_names=[output_var_name],
        input_levels=[input_level] * num_inputs,
        output_levels=[output_level] * num_inputs
    )

    # Script the wrapper (not trace) to preserve attributes
    scripted_model = torch.jit.script(wrapper)

    # Save the model
    scripted_model.save(output_path)

    print(f"Created emulator: {output_path}")
    print(f"  Output variable: {output_var_name}")
    print(f"  Input variables: {input_var_names}")
    print(f"  Jacobian values: {jacobian_values}")
    print(f"  Input level: {input_level}, Output level: {output_level}")

    return scripted_model


class VerticalLinearEmulator(nn.Module):
    """
    Vertical emulator for testing TorchBalance vertical emulators.

    Uses a constant diagonal Jacobian: J[k_out, k_in] = diag_value if k_out == k_in else 0.
    jac_physical returns [nnodes, nRequestedPairs] as required by setupVerticalEmulator.
    No input_levels / output_levels attributes — the C++ reads level counts from Atlas field shapes.
    """

    def __init__(self, n_out_levels: int, n_in_levels: int, diag_value: float = 0.5):
        super().__init__()
        jac = torch.zeros(n_out_levels, n_in_levels)
        for k in range(min(n_out_levels, n_in_levels)):
            jac[k, k] = diag_value
        self.register_buffer('jac_matrix', jac)
        self.n_out_levels = n_out_levels
        self.n_in_levels = n_in_levels

    def forward(self, inputs: torch.Tensor) -> torch.Tensor:
        return torch.matmul(inputs, self.jac_matrix.t())

    @torch.jit.export
    def jac_physical(self, inputs: torch.Tensor, mask: torch.Tensor,
                     requested_row_indices: torch.Tensor,
                     requested_col_indices: torch.Tensor) -> torch.Tensor:
        """
        Returns only requested Jacobian entries: [nnodes, nRequestedPairs].

        requested_row_indices: 1-D int64 tensor of output level indices.
        requested_col_indices: 1-D int64 tensor of input level indices into the full input
        level dimension.
        """
        nnodes = inputs.shape[0]
        req_jac = self.jac_matrix[requested_row_indices, requested_col_indices]
        jac = req_jac.unsqueeze(0).expand(nnodes, -1).clone()
        return jac * mask


def create_vertical_emulator_with_metadata(
    output_var_name,
    input_var_names,
    n_out_levels,
    n_in_levels,
    diag_value=0.5,
    output_path="vertical_emulator.ts"
):
    """
    Create a TorchScript vertical emulator for TorchBalance vertical emulator tests.

    The emulator has a constant diagonal Jacobian. No input_levels / output_levels
    attributes are stored — setupVerticalEmulator derives level counts from Atlas fields.
    """
    model = VerticalLinearEmulator(n_out_levels=n_out_levels,
                                   n_in_levels=n_in_levels,
                                   diag_value=diag_value)
    model.eval()

    class VerticalEmulatorWithMetadata(torch.nn.Module):
        def __init__(self, inner, input_names, output_names):
            super().__init__()
            self.model = inner
            self.input_names = torch.jit.Attribute(input_names, List[str])
            self.output_names = torch.jit.Attribute(output_names, List[str])

        def forward(self, x: torch.Tensor) -> torch.Tensor:
            return self.model(x)

        @torch.jit.export
        def jac_physical(self, x: torch.Tensor, mask: torch.Tensor,
                         requested_row_indices: torch.Tensor,
                         requested_col_indices: torch.Tensor) -> torch.Tensor:
            return self.model.jac_physical(x, mask, requested_row_indices, requested_col_indices)

    wrapper = VerticalEmulatorWithMetadata(
        model,
        input_names=input_var_names,
        output_names=[output_var_name],
    )
    scripted = torch.jit.script(wrapper)
    scripted.save(output_path)

    print(f"Created vertical emulator: {output_path}")
    print(f"  Output variable: {output_var_name} ({n_out_levels} levels)")
    print(f"  Input variables: {input_var_names} ({n_in_levels} levels each)")
    print(f"  Jacobian: diagonal {diag_value}")

    return scripted


def create_test_suite():
    """
    Create a complete test suite with multiple emulators.

    Creates simple linear emulators for testing TorchBalance infrastructure.
    The models are intentionally simple for predictable, verifiable behavior.

    Creates two emulators:
    1. air_temperature depends on [sst, qv]
    2. sea_water_potential_temperature depends on [tair, qv]
    3. water_vapor_mixing_ratio_wrt_moist_air has no emulator (identity only)
    """
    print("Creating test suite of emulators...\n")

    # Emulator 1: air_temperature
    create_emulator_with_metadata(
        output_var_name="air_temperature",
        input_var_names=["sea_water_potential_temperature", "water_vapor_mixing_ratio_wrt_moist_air"],
        jacobian_values=[0.3, 0.2],  # dT/dSST = 0.3, dT/dQV = 0.2
        input_level=0,
        output_level=0,
        output_path="air_temperature_emulator.ts"
    )
    print()

    # Emulator 2: sea_water_potential_temperature
    create_emulator_with_metadata(
        output_var_name="sea_water_potential_temperature",
        input_var_names=["air_temperature", "water_vapor_mixing_ratio_wrt_moist_air"],
        jacobian_values=[0.4, 0.1],  # dSST/dT = 0.4, dSST/dQV = 0.1
        input_level=0,
        output_level=0,
        output_path="sea_water_potential_temperature_emulator.ts"
    )
    print()

    print("\nTest suite created!")
    print("\nExpected matrix structure:")
    print("                       sst    qv    tair")
    print("  sst                [  I    0.1   0.4  ]")
    print("  qv                 [  0     I     0   ]")
    print("  tair               [ 0.3   0.2    I   ]")
    print("\nYAML configuration:")
    print("  surface emulators:")
    print("    - name: air_temperature")
    print("      path: air_temperature_emulator.ts")
    print("      jacobian wrt: [sea_water_potential_temperature, water_vapor_mixing_ratio_wrt_moist_air]")
    print("    - name: sea_water_potential_temperature")
    print("      path: sea_water_potential_temperature_emulator.ts")
    print("      jacobian wrt: [air_temperature, water_vapor_mixing_ratio_wrt_moist_air]")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create TorchScript emulators for testing TorchBalance")
    parser.add_argument("output", nargs="?", default="test_emulator.ts",
                        help="Output path for test emulator (default: test_emulator.ts)")
    args = parser.parse_args()

    # Derive output directory from the provided path
    output_dir = Path(args.output)

    # Create a complete test suite with multiple emulators
    print("Creating test suite of emulators...\n")

    # Emulator 1: air_temperature
    create_emulator_with_metadata(
        output_var_name="air_temperature",
        input_var_names=[
            "sea_water_potential_temperature",
            "water_vapor_mixing_ratio_wrt_moist_air"
        ],
        jacobian_values=[0.3, 0.2],
        input_level=0,
        output_level=0,
        output_path=str(
            output_dir / "air_temperature_emulator.ts"
        )
    )
    print()

    # Emulator 2: sea_water_potential_temperature
    create_emulator_with_metadata(
        output_var_name="sea_water_potential_temperature",
        input_var_names=[
            "air_temperature",
            "water_vapor_mixing_ratio_wrt_moist_air"
        ],
        jacobian_values=[0.4, 0.1],
        input_level=0,
        output_level=0,
        output_path=str(
            output_dir / "sst_emulator.ts"
        )
    )
    print()

    # Vertical emulator: sea_water_salinity from sea_water_potential_temperature
    # Diagonal Jacobian 0.5 over 4 levels (matches CS-LFR-15 / 4-level test geometry)
    create_vertical_emulator_with_metadata(
        output_var_name="sea_water_salinity",
        input_var_names=["sea_water_potential_temperature"],
        n_out_levels=4,
        n_in_levels=4,
        diag_value=0.5,
        output_path=str(output_dir / "vertical_emulator.ts")
    )
    print()

    print("\nTest suite created!")
    print("\nExpected matrix structure (surface):")
    print("                       sst    qv    tair")
    print("  sst                [  I    0.1   0.4  ]")
    print("  qv                 [  0     I     0   ]")
    print("  tair               [ 0.3   0.2    I   ]")
    print()
    print("\nVertical emulator (sea_water_salinity w.r.t. sea_water_potential_temperature):")
    print("  J_sal_temp = 0.5 * I_4  (4x4 diagonal Jacobian block)")
    print()
    print("  Full K operator is 8x8, partitioned by variable (temp: levels 0-3, sal: levels 0-3):")
    print()
    print("              temp (4 lvls)   sal (4 lvls)")
    print("  temp (4 lvls) [    I_4           0     ]")
    print("  sal (4 lvls) [  0.5*I_4         I_4   ]")
    print()
