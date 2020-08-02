//
//  Renderer.swift
//  MetalBench
//
//  Created by Alia on 30/07/2020.
//

import AppKit
import MetalKit

enum Scene: Int, CaseIterable {
	case glass = 0, forest = 1
	
	var name: String {
		get {
			switch self {
				case .glass:
					return "Glass"
				case .forest:
					return "Forest floor"
			}
		}
	}
	
	var isRealtime: Bool {
		get {
			switch self {
				case .glass:
					return true
				case .forest:
					return false
			}
		}
	}
}

class Renderer: NSObject, ObservableObject {
	@Published var gpuList: [MTLDevice]
	@Published var selectedGPU = 0 {
		didSet {
			dev = gpuList[selectedGPU]
		}
	}
	
	@Published var sceneList: [Scene]
	@Published var selectedScene = 0 {
		didSet {
			isRealtime = sceneList[selectedScene].isRealtime
		}
	}
	var isRealtime = true {
		didSet {
			currentRayCount = isRealtime ? 1 : 2
		}
	}
	
	@Published var fps = 0.0
	@Published var megaRaysPerSecond = 0
	@Published var averageMegaRaysPerSecond = 0
	
	var last30RayCounts = Array(repeating: 0, count: 30)
	var count30 = 0 // Incremented when adding a ray count value, average is appended to previous averages when counter hits 30
	var previousRayCountAverages = [Int]()
	var wasReset = true // When this is true, last30RayCounts is populated with the next value
	
	var currentRayCount = 1
	
	var dev: MTLDevice
	
	override init() {
		print("Initialising preview")
		
		// Get GPU list, use first GPU
		let gpus = MTLCopyAllDevices()
		guard gpus.count > 0 else {
			exit(EXIT_FAILURE)
		}
		
		gpuList = gpus
		dev = gpus.first!
		
		sceneList = Scene.allCases
	}
	
	func addFrameDuration(_ duration: CFTimeInterval) {
		guard duration > 0.0, currentRayCount > 0 else {
			return
		}
		
		// If we're not realtime then duration represents 1/32 of a frame, scale accordingly
		let totalDuration = duration * (isRealtime ? 1.0 : 32.0)
		
		// FPS is live value
		fps = 1.0 / totalDuration
		
		// RPS is based on resolution * rays divided by time taken
		let raysPerSecond = Int(CFTimeInterval(1280 * 720 * currentRayCount) / totalDuration)
		
		// Add to last 30
		if wasReset {
			// After reset, the short term history is filled with the initial value
			last30RayCounts = Array(repeating: raysPerSecond, count: 30)
			wasReset = false
		} else {
			last30RayCounts[count30] = raysPerSecond
		}
		count30 += 1
		
		// If last 30 is full, average and append to long term history
		if count30 == 30 {
			let sum = last30RayCounts.reduce(0, +)
			previousRayCountAverages.append(sum / 30)
			count30 = 0
		}
		
		// Find short term average ray count
		let currentAverage = last30RayCounts.reduce(0, +) / 30
		
		// The standard MRAYS value is based on short term average
		megaRaysPerSecond = currentAverage / 1_000_000
		
		// Average rays is average over complete history
		// Add current average to sum of previous averages, divide by count
		averageMegaRaysPerSecond = ((previousRayCountAverages.reduce(0, +) + currentAverage) / (previousRayCountAverages.count + 1)) / 1_000_000
		
		// If realtime, calculate optimum ray count for target fps
		if sceneList[selectedScene].isRealtime {
			let raysPerPixel = Int(floor(Double(currentAverage) / (1280.0 * 720.0 * 30.0)))
			currentRayCount = max(1, raysPerPixel)
		}
	}
	
	func resetStats() {
		count30 = 0
		previousRayCountAverages.removeAll()
		wasReset = true
	}
}
