//
//  PreviewView.swift
//  MetalBench
//‪
//  Created by Alia on 30/07/2020.
//

import SwiftUI
import MetalKit

struct PreviewView: NSViewRepresentable {
	typealias NSViewType = MTKView
	
	@EnvironmentObject var renderer: Renderer
	
	func makeCoordinator() -> PreviewCoordinator {
		PreviewCoordinator(self)
	}
	
	func makeNSView(context: NSViewRepresentableContext<PreviewView>) -> MTKView {
		let mtkView = MTKView()
		mtkView.delegate = context.coordinator
		mtkView.device = renderer.dev
		mtkView.autoResizeDrawable = false
		mtkView.drawableSize = CGSize(width: 1280, height: 720)
		mtkView.framebufferOnly = false
		
		return mtkView
	}
	
	func updateNSView(_ nsView: MTKView, context: NSViewRepresentableContext<PreviewView>) {
		// Check view is set
		if context.coordinator.view == nil {
			context.coordinator.view = nsView
		}
		
		// Check for device change
		let activeDev = context.coordinator.parent.renderer.dev
		if activeDev.registryID != context.coordinator.metalDevice.registryID {
			// New device, update the view and coordinator
			context.coordinator.metalDevice = activeDev
			nsView.device = activeDev
			context.coordinator.configure()
			context.coordinator.parent.renderer.resetStats()
		}
		
		// Check for scene change
		let activeScene = context.coordinator.parent.renderer.selectedScene
		if activeScene != context.coordinator.currentScene {
			context.coordinator.configure()
			context.coordinator.parent.renderer.resetStats()
		}
	}
	
}

class PreviewCoordinator : NSObject, MTKViewDelegate {
	var parent: PreviewView
	var view: MTKView?
	var metalDevice: MTLDevice
	var queue: MTLCommandQueue!
	var computePS: MTLComputePipelineState?
	var displayPS: MTLComputePipelineState?
	
	var currentScene = 0 {
		didSet {
			isRealtime = parent.renderer.sceneList[currentScene].isRealtime
		}
	}
	var isRealtime = true {
		didSet {
			// Either display-driven rendering or frame-available driven depending on scene
			view?.isPaused = !isRealtime
			view?.enableSetNeedsDisplay = !isRealtime
			
			if !isRealtime {
				// Create the storage texture for incremental rendering
				let desciptor = MTLTextureDescriptor.texture2DDescriptor(pixelFormat: .rgba32Float, width: 1280, height: 720, mipmapped: false)
				desciptor.usage = [.shaderRead, .shaderWrite]
				desciptor.storageMode = .private
				incrementalTexture = metalDevice.makeTexture(descriptor: desciptor)
			}
		}
	}
	var endIncrementalRendering = false
	var incrementalTexture: MTLTexture?
	
	// Incremental renders tiled, 32x2 tiles
	var incrementalTile = SIMD2<UInt>(0, 0)
	
	var viewConfigured = false
	
	var rayCountBuffer: MTLBuffer?
	var gridOffsetBuffer: MTLBuffer?
	var timeBuffer: MTLBuffer?
	var startTime = Date()
	
	let targetFPS = 30.0
	
	init(_ parent: PreviewView) {
		self.parent = parent
		self.metalDevice = parent.renderer.dev
		super.init()
	}
	
	/// Prepares rendering
	func configure() {
		endIncrementalRendering = true
		guard let library = metalDevice.makeDefaultLibrary() else {
			print("Failed to create libary")
			return
		}
		
		let constants = MTLFunctionConstantValues()
		var value = Int32(parent.renderer.selectedScene)
		constants.setConstantValue(&value, type: MTLDataType.int, index: 0)
		
		let kernel: MTLFunction
		do {
			kernel = try library.makeFunction(name: "k", constantValues: constants)
		} catch let e {
			print("Failed to create kernel")
			print(e.localizedDescription)
			return
		}
		
		do {
			computePS = try metalDevice.makeComputePipelineState(function: kernel)
		} catch let e {
			print ("Failed to make compute pipeline")
			print(e.localizedDescription)
			return
		}
		
		currentScene = parent.renderer.selectedScene
		
		queue = metalDevice.makeCommandQueue()
		
		// Set up the buffers to contain time and resolution
		startTime = Date()
		timeBuffer = metalDevice.makeBuffer(length: MemoryLayout<Float>.size, options: [])
		rayCountBuffer = metalDevice.makeBuffer(length: MemoryLayout<Int>.size, options: [])
		gridOffsetBuffer = metalDevice.makeBuffer(length: MemoryLayout<SIMD2<UInt>>.size, options: [])
		updateRayCount()
		
		viewConfigured = true
		
		if !isRealtime {
			// Create a display pipeline
			let displayKernel: MTLFunction
			do {
				displayKernel = try library.makeFunction(name: "displayIntegrated", constantValues: constants)
			} catch let e {
				print("Failed to create kernel")
				print(e.localizedDescription)
				return
			}
			
			do {
				displayPS = try metalDevice.makeComputePipelineState(function: displayKernel)
			} catch {
				print ("Failed to make display pipeline")
				return
			}
			endIncrementalRendering = false
			incrementalTile = SIMD2<UInt>(0, 0)
			drawIncremental()
		}
	}
	
	func updateTime() {
		guard let buffer = self.timeBuffer else { return }
		
		var t = Float(Date().timeIntervalSince(startTime) as Double)
		
		let bufferPointer = buffer.contents()
		memcpy(bufferPointer, &t, MemoryLayout<Float>.size)
	}
	
	func updateRayCount() {
		guard let buffer = self.rayCountBuffer else { return }
		
		let bufferPointer = buffer.contents()
		var rays = isRealtime ? parent.renderer.currentRayCount : 1
		memcpy(bufferPointer, &rays, MemoryLayout<Int>.size)
	}
	
	func updateGridOffset() {
		guard let buffer = self.gridOffsetBuffer else { return }
		
		let bufferPointer = buffer.contents()
		memcpy(bufferPointer, &incrementalTile, MemoryLayout<Int>.size)
	}
	
	func mtkView(_ view: MTKView, drawableSizeWillChange size: CGSize) {
	}
	
	/// Renders incrementally into a backing texture, and marks the view as needing an update
	func drawIncremental() {
		updateTime()
		updateGridOffset()
		
		let commandBuffer = queue.makeCommandBuffer()
		
		guard let pipe = computePS, let buffer = commandBuffer, let cmdEncoder = buffer.makeComputeCommandEncoder() else { return }
		
		buffer.addCompletedHandler { cmdBuffer in
			// Calculate GPU side frame duration
			let start = cmdBuffer.gpuStartTime
			let end = cmdBuffer.gpuEndTime
			let duration = end - start
			
			// Update state for UI update
			DispatchQueue.main.async {
				self.parent.renderer.addFrameDuration(duration)
				self.view?.setNeedsDisplay(NSRect(origin: .zero, size: CGSize(width: 1280, height: 720)))
			}
			
			// Draw again if not ended
			if !self.endIncrementalRendering { self.drawIncremental() }
		}
		#warning("Crash on change GPU in forest due to pipe associated with different GPU")
		cmdEncoder.setComputePipelineState(pipe)
		
		cmdEncoder.setTexture(incrementalTexture, index: 0)
		cmdEncoder.setTexture(incrementalTexture, index: 1)
		cmdEncoder.setBuffer(timeBuffer!, offset: 0, index: 0)
		cmdEncoder.setBuffer(rayCountBuffer!, offset: 0, index: 1)
		cmdEncoder.setBuffer(gridOffsetBuffer!, offset: 0, index: 2)
		
		let threadGroupCount = MTLSizeMake(8, 8, 1)
		let threadGroups = MTLSizeMake((1280 / 32) / threadGroupCount.width, (720 / 1) / threadGroupCount.height, 1)
		cmdEncoder.dispatchThreadgroups(threadGroups, threadsPerThreadgroup: threadGroupCount)
		cmdEncoder.endEncoding()
		
		buffer.commit()
		
		incrementalTile.x += 1280 / 32
		if incrementalTile.x >= 1280 { incrementalTile.x = 0 }
	}
	
	func draw(in view: MTKView) {
		if !viewConfigured { configure() }
		
		guard let drawable = view.currentDrawable else {
			return
		}
		
		updateTime()
		updateRayCount()
		
		let commandBuffer = queue.makeCommandBuffer()
		
		guard let pipe = isRealtime ? computePS : displayPS, let buffer = commandBuffer, let cmdEncoder = buffer.makeComputeCommandEncoder() else { return }
		
		if isRealtime {
			buffer.addCompletedHandler { cmdBuffer in
				// Calculate GPU side frame duration
				let start = cmdBuffer.gpuStartTime
				let end = cmdBuffer.gpuEndTime
				let duration = end - start
				
				// Update state for UI update
				DispatchQueue.main.async {
					self.parent.renderer.addFrameDuration(duration)
				}
			}
		}
		
		cmdEncoder.setComputePipelineState(pipe)
		
		if isRealtime {
			cmdEncoder.setTexture(drawable.texture, index: 0)
			cmdEncoder.setBuffer(timeBuffer!, offset: 0, index: 0)
			cmdEncoder.setBuffer(rayCountBuffer!, offset: 0, index: 1)
		} else {
			cmdEncoder.setTexture(incrementalTexture, index: 0)
			cmdEncoder.setTexture(drawable.texture, index: 1)
		}
		
		let threadGroupCount = MTLSizeMake(8, 8, 1)
		let threadGroups = MTLSizeMake(drawable.texture.width / threadGroupCount.width, drawable.texture.height / threadGroupCount.height, 1)
		cmdEncoder.dispatchThreadgroups(threadGroups, threadsPerThreadgroup: threadGroupCount)
		cmdEncoder.endEncoding()
		
		buffer.present(drawable)
		buffer.commit()
	}
}
