//
//  Modifiers.swift
//  MetalBench
//
//  Created by Alia on 30/07/2020.
//

import SwiftUI

struct UITextModifier: ViewModifier {
	func body(content: Content) -> some View {
		content
			.foregroundColor(Color.white)
	}
}
