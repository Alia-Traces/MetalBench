//
//  ContentView.swift
//  MetalBench
//
//  Created by Alia on 30/07/2020.
// 

import SwiftUI

struct ContentView: View {
	@EnvironmentObject var renderer: Renderer
	
	var body: some View {
		ZStack(alignment: .bottom) {
			PreviewView()
				.frame(width: 1280, height: 720, alignment: .center)
			StatusOverlay()
		}
		.frame(width: 1280, height: 720, alignment: .center)
	}
}

struct ContentView_Previews: PreviewProvider {
	static let renderer = Renderer()
	
	static var previews: some View {
		ContentView()
			.environmentObject(renderer)
	}
}
