/*******************************************************************************
 * Copyright (c) 2015 David Lamparter, Daniel Marbach
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *******************************************************************************/
package ch.unil.genescore.vegas;

import com.sun.jna.Library;
import com.sun.jna.ptr.DoubleByReference;
import com.sun.jna.ptr.IntByReference;

public interface MvnPack extends Library {

		/**
		 * See http://www.math.wsu.edu/faculty/genz/software/fort77/mvtdstpack.f
		 * @param n
		 * @param lower
		 * @param upper
		 * @param infin
		 * @param correl
		 * @param maxpts
		 * @param abseps
		 * @param releps
		 * @param error
		 * @param value
		 * @param inform
		 */
		void mvtdst_(IntByReference n,IntByReference df, double[] lower, double[] upper, int[] infin, double[] correl, double[] delta,
				IntByReference maxpts, DoubleByReference abseps, DoubleByReference releps, 
				DoubleByReference error, DoubleByReference value, IntByReference inform);
	
		/**
		 * See http://www.math.wsu.edu/faculty/genz/software/fort77/mvnexppack.f
		 * @param n
		 * @param lower
		 * @param upper
		 * @param infin
		 * @param correl
		 * @param maxpts
		 * @param abseps
		 * @param releps
		 * @param error
		 * @param value
		 * @param inform
		 */
	//	void mvnexp_(IntByReference n, double[] lower, double[] upper, int[] infin, double[] correl,
		//		IntByReference maxpts, DoubleByReference abseps, DoubleByReference releps,
			//	double[] error, double[] value, IntByReference inform);

		/**
		 * See http://www.math.wsu.edu/faculty/genz/software/fort77/mvnxpppack.f
		 * @param n
		 * @param lower
		 * @param upper
		 * @param infin
		 * @param correl
		 * @param maxpts
		 * @param abseps
		 * @param releps
		 * @param error
		 * @param value
		 * @param inform
		 */
	//	void mvnxpp_(IntByReference n, double[] lower, double[] upper, int[] infin, double[] correl,
		//		IntByReference maxpts, DoubleByReference abseps, DoubleByReference releps,
			//	double[] error, double[] value, IntByReference inform);

}


