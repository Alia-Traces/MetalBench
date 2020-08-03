//
//  MtlDeviceExtensions.swift
//  MetalBench
//
//  Created by Alia on 30/07/2020.
//

import Metal

extension MTLDevice {
	var uiDescription: String {
		get {
			return self.name + (self.isLowPower ? " (low power)" : "") + (self.isRemovable ? " (eGPU)" : "")
		}
	}
}
