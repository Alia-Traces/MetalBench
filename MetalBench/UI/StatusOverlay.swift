//
//  StatusOverlay.swift
//  MetalBench
//
//  Created by Alia on 30/07/2020.
//

import SwiftUI

struct StatusOverlay: View {
	@EnvironmentObject var renderer: Renderer
	
	var body: some View {
		ZStack {
			VStack {
				HStack {
					Picker("GPU", selection: $renderer.selectedGPU) {
						ForEach(0..<renderer.gpuList.count) { index in
                            Text(self.renderer.gpuList[index].uiDescription).tag(index)
						}
					}
					.frame(width: 400, height: nil, alignment: .center)
					
					Spacer()
					
					Text(String(format: "%1.0f", renderer.fps))
					Text("fps")
					
				}
				HStack {
					Picker("Scene", selection: $renderer.selectedScene) {
						ForEach(0..<renderer.sceneList.count) { index in
                            Text(self.renderer.sceneList[index].name).tag(index)
						}
					}
					.frame(width: 400, height: nil, alignment: .center)
					
					Spacer()
					
					Text("\(renderer.megaRaysPerSecond)")
					Text("MRays / Second")
				}
				HStack {
					Text("1280 x 720")
					
					Spacer()
					
					Text("\(renderer.averageMegaRaysPerSecond)")
					Text("MRays / Second average")
				}
			}
			.padding(10)
			.background(Color(.sRGB, white: 0.0, opacity: 0.25))
		}
		.modifier(UITextModifier())
	}
}

struct StatusOverlay_Previews: PreviewProvider {
	static let renderer = Renderer()
	
	static var previews: some View {
		StatusOverlay()
			.frame(width: 1280, height: nil, alignment: .center)
			.environmentObject(renderer)
	}
}
